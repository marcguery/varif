from sys import stderr
from .config import Config
from .vcf import Vcf, Vcfchunk
from .variants import Variants
from .families import Families
from .annotations import Annotations
from .fasta import Fasta
from .version import __version__
import time
import math
from multiprocessing import Pool


class Connection(object):
    """Filter and annotate variants with samples optionnaly grouped by 
    familiy or lineage. Will do as many passes as there are combinations to compare
    
    """

    def __init__(self):
        """
        vcf (Vcf) : Vcf data
        families (Families) : PED data
        annotations (Annotations) : GFF data and tools
        fasta (Fasta) : Fasta object
        samples (list) : Sorted list of samples processed
        group1 (list) : Samples belonging to the 'group1' to be compared
        group2 (list) : Samples belonging to the 'group2' to be compared
        chunks (int) : Number of chunks to use to process the variants separately
        comparison (str) : Type of comparison to use between groups
        intergenicRegions (bool) : Whether to include intergenic regions in the output
        outFile (str) : Name of the output file name
        outputVcf (bool) : Whether to also output the filtered VCF file

        """

        start_time = time.time()
        Config.set_options()
        if Config.options["version"]:
            print(__version__)
            raise SystemExit
        Config.check_options()
        print(("Running Varif %s"
              " with options \n%s")%(__version__,
                                     "\n".join(" : ".join([Config.long_options[opt], str(Config.options[opt])]) for opt in Config.options)))
        self.vcf = None
        self.families = Families()
        if Config.options["ped"] is not None:
            self.families.read_ped(Config.options["ped"])
        self.annotations = Annotations()
        self.annotations.load_annotations_from_GFF(Config.options["gff"])
        self.fasta = Fasta()
        self.fasta.load_data_from_FASTA(Config.options["fasta"])
        Config.verbose_print("Genome sequences, annotation and sample metadata loaded in %s seconds"%(round(time.time() - start_time)))

        self.samples = []
        self.group1 = []
        self.group2 = []
        self.chunks = 1
        self.comparison = Config.options["comparison"]
        self.intergenicRegions = Config.options["intergenicRegions"]
        self.outFile = Config.options["outFile"]
        self.outputVcf = Config.options["outputVcf"]
        self.get_all_groups()
        print("Varif run was completed in %s seconds"%(round(time.time() - start_time)))
        
    
        
    def print_log(self, allvariants, variantId):
        """
        Print the log stored in a variant if not empty

        allvariants (Variants) : Current variants processed
        variantId (str) : The unique identifier of the variant

        """
        log = ""
        mainlog = allvariants.variants[variantId]["log"][0]
        varlogs = allvariants.variants[variantId]["log"][1]

        if mainlog != "":
            log = "Variant %s: %s"%(variantId, mainlog)
        if "".join(varlogs) != "": 
            pass
        if log != "":
            print(log, file = stderr)
    
    def get_aa_from_mutation(self, chromosome, position, reference, mutation, gffId, stripMutation = False):
        """
        Get the change of aminoacid induced by a mutation

        chromosome (str) : Chromosome where the mutation occurs
        position (int) : Position of the mutation
        reference (str) : Upper case sequence before mutation
        mutation (str) : Upper case sequence after mutation
        gffId (str) : ID of the CDS as in GFF3 file
        stripMutation (bool) : Whether to remove bases of the mutated sequence located in introns

        return (list) : CDS and aminoacids (with their positions) before and after mutation

        """
        assert self.annotations.annotations[gffId]['annotation'] == 'CDS'
        strand = self.annotations.annotations[gffId]['strand']
        windowBefore = 2+Config.options["protWindowBefore"]*3
        windowAfter = 2+Config.options["protWindowAfter"]*3

        if len(self.annotations.annotations[gffId]["parents"]) > 1:
            print("Ambiguous origin of CDS %s: one of %s. Taking the first one (%s) by default"%(gffId, ", ".join(self.annotations.annotations[gffId]["parents"]), self.annotations.annotations[gffId]["parents"][0]), file = stderr)
        allCDSsorted, cdsCoords = self.annotations.get_CDS_from_same_parent(gffId)
        intronCoords = self.fasta.get_introns_coords(cdsCoords)
        genesequence = self.fasta.merge_CDS(chromosome, cdsCoords)
        
        if len(reference) > 1:
            #Remove bases from the reference if located in an intronic region
            newReference, newReferencePosition = self.fasta.remove_bases_from_features(reference, intronCoords, position-1)
            assert newReference != "", "The reference sequence is not located in a CDS"
        else:
            #For SNPs, no difference : if they were located in an intron, this function would not be called
            newReference, newReferencePosition = [reference, position - 1]
        
        #Remove bases from the mutated sequence only if stripMutation is enabled :
        # the resulting CDS could be missing actual codons (coming from a genuine CDS region)
        # but false positive will be removed (codons coming from an intronic region)
        if stripMutation:
            newMutation, newMutationPosition = self.fasta.remove_bases_from_features(mutation, intronCoords, position-1)
        else:
            newMutation = mutation
        
        #Relative position of the reference sequence to the start position of the gene
        relativeVariantPosition = self.fasta.get_relative_gene_position(cdsCoords, newReferencePosition)
        startIndex = max(0, relativeVariantPosition)
        endIndex = min(len(genesequence), startIndex+len(newReference))
        oldCDS, aaPos, phase = self.fasta.window_CDS(genesequence, strand, startIndex, endIndex, windowBefore, windowAfter)
        oldProt = self.fasta.translate_CDS(oldCDS, strand, phase)

        if newMutation == newReference:
            #If both original and mutated sequences (without intron bases) are identical,
            # no need to translate twice the same sequence
            aaChanges = [[aaPos,int(len(genesequence)/3)], oldCDS, oldProt, [0,0], "", ""]    
        else:
            genesequencemut, startIndexMut, endIndexMut = self.fasta.insert_mutation(chromosome, genesequence, strand, cdsCoords[0][0], cdsCoords[-1][1], relativeVariantPosition, newReference, newMutation)
            newCDS, newaaPos, newphase = self.fasta.window_CDS(genesequencemut, strand, startIndexMut, endIndexMut, windowBefore, windowAfter)
            newProt = self.fasta.translate_CDS(newCDS, strand, newphase)
            aaChanges = [[aaPos,int(len(genesequence)/3)], oldCDS, oldProt, [newaaPos,int(len(genesequencemut)/3)], newCDS, newProt]   
        
        return aaChanges
    
    def load_features(self, allvariants, variantId):
        """
        Store aminoacids obtained from CDS translations
        
        allvariants (Variants) : Current variants processed
        variantId (str) : Unique identifier of the variant

        """
        for gffId in allvariants.variants[variantId]["features"]:
            if self.annotations.annotations[gffId]['annotation'] != 'CDS':
                continue #Not a CDS feature
            allvariants.variants[variantId]["aaPosAlts"][gffId] = []
            allvariants.variants[variantId]["aaAlts"][gffId] = []
            allvariants.variants[variantId]["cdsAlts"][gffId] = []
            for alt in allvariants.variants[variantId]["alts"]:
                aaDiff = self.get_aa_from_mutation(allvariants.variants[variantId]["chromosome"], 
                allvariants.variants[variantId]["position"], 
                allvariants.variants[variantId]["ref"], alt, gffId)
                allvariants.variants[variantId]["aaPosAlts"][gffId].append(aaDiff[3])
                allvariants.variants[variantId]["cdsAlts"][gffId].append(aaDiff[4])
                allvariants.variants[variantId]["aaAlts"][gffId].append(aaDiff[5])
            allvariants.variants[variantId]["aaPosRef"][gffId] = aaDiff[0]
            allvariants.variants[variantId]["cdsRef"][gffId] = aaDiff[1]
            allvariants.variants[variantId]["aaRef"][gffId] = aaDiff[2]
    
    def define_categories(self, allvariants, variantId):
        """
        Determine if the variant alleles are a SNP or INDEL and add genomic location

        allvariants (Variants) : Current variants processed
        variantId (str) : Unique identifier of the variant

        """
        ref = allvariants.variants[variantId]["ref"]
        for alt in allvariants.variants[variantId]["alts"]:
            if len(ref) == len(alt):
                if ref[0] != alt[0] and ref[1:] == alt[1:]:
                    maincategory = "SNP"
                else:
                    raise ValueError("VCF SNP is not formatted properly, unexpected REF (%s) and ALT (%s)")
            else:
                maincategory = "INDEL"
            subcategories = sorted(set([self.annotations.annotations[gffId]['annotation'] for gffId in allvariants.variants[variantId]["features"]]))
            allvariants.variants[variantId]["categories"].append(",".join([maincategory]+[subcategories for subcategories in subcategories]))

    def print_header(self, csvsep = ";"):
        """
        Print the header of the output file(s)

        csvsep (str) : Separator to use for the CSV file

        return (list) : The headers of each of the CSV and filtered VCF files

        """
        baseHeader = [
                "Chromosome", "Position", "Type", "Ref", "Alt", 
                "CDS_ref", "CDS_alt", "AA_ref", "AA_alt", "Annotation", 
                "Missingness_G1", "Missingness_G2",
                "Heterozygous_G1", "Heterozygous_G2",
                "Alt_frequency_G1", "Alt_frequency_G2"]
        csvline = csvsep.join(baseHeader+self.samples)+"\n"

        vcfline = "".join(self.vcf.vcfheader)
        return [csvline, vcfline]
    
    def annotate_variant(self, allvariants, variantId = None, altIndex = 0, sep = ";"):
        """
        Annotate a single VCF variant for the main output of varif separating each alt's variant

        allvariants (Variants) : Current variants processed
        variantId (str) : The unique identifier of the variant
        altIndex (int) : The index of the alt of this variant
        sep (str) : The field separator

        return (str) : The line of the header or of each alt

        """
        assert variantId is not None

        windowLenSum = len(allvariants.variants[variantId]["refwindow"][0])+len(allvariants.variants[variantId]["refwindow"][1])
        leftwindow = allvariants.variants[variantId]["refwindow"][0]+"|" if windowLenSum > 0 else ""
        rightwindow = "|"+allvariants.variants[variantId]["refwindow"][1] if windowLenSum > 0 else ""
        content = [
            allvariants.variants[variantId]["chromosome"],                              #0                     
            str(allvariants.variants[variantId]["position"]),                           #1
            allvariants.variants[variantId]["categories"][altIndex],                    #2
            leftwindow+allvariants.variants[variantId]["ref"]+rightwindow,              #3
            allvariants.variants[variantId]["alts"][altIndex],                          #4
            "NA", "NA", "NA","NA", "NA",                                                #5-9
            '{:.6f}'.format(allvariants.variants[variantId]["miss"][altIndex][0]),      #10
            '{:.6f}'.format(allvariants.variants[variantId]["miss"][altIndex][1]),      #11
            '{:.6f}'.format(allvariants.variants[variantId]["heteroz"][altIndex][0]),   #12
            '{:.6f}'.format(allvariants.variants[variantId]["heteroz"][altIndex][1]),   #13
            '{:.6f}'.format(allvariants.variants[variantId]["props"][altIndex][0]),     #14
            '{:.6f}'.format(allvariants.variants[variantId]["props"][altIndex][2])      #15
            ]
        content += ["NA"]*len(allvariants.samples)

        aaposref = []
        cdsref = []
        aaref = []
        aaposalts = []
        cdsalts = []
        aaalts = []
        for gffId in allvariants.variants[variantId]["cdsRef"]:
            aaposref.append("/".join(str(pos) for pos in allvariants.variants[variantId]["aaPosRef"][gffId]))
            cdsref.append(allvariants.variants[variantId]["cdsRef"][gffId]+" ("+gffId+")")
            aaref.append(allvariants.variants[variantId]["aaRef"][gffId])
            aaposalts.append("/".join(str(pos) for pos in allvariants.variants[variantId]["aaPosAlts"][gffId][altIndex]))
            cdsalts.append(allvariants.variants[variantId]["cdsAlts"][gffId][altIndex])
            aaalts.append(allvariants.variants[variantId]["aaAlts"][gffId][altIndex])
            
        if len(allvariants.variants[variantId]["cdsRef"]) > 0: #There is a feature at least
            content[5] = ":".join(cdsref) #the same cds ref in different features
            content[6] = ":".join(cdsalts) #different cds alt
            content[7] = "):".join(" (".join(str(cont) for cont in couple) for couple in zip(aaref, aaposref))+")" #the same aa ref in different features
            content[8] = "):".join(" (".join(str(cont) for cont in couple) for couple in zip(aaalts, aaposalts))+")" #different aa alt
        annotations = set(self.annotations.annotations[geneId]['description']+" ("+geneId+")" for geneId in allvariants.variants[variantId]["features"] if "gene" in self.annotations.annotations[geneId]['annotation'])
        content[9] = ":".join(annotations) if annotations != set() else content[9] #potentially different annotations
        for i, sample in enumerate(allvariants.samples):#asps
            content[16+i] = '{:.6f}'.format(allvariants.variants[variantId]["asps"][sample][altIndex+1])
        line = sep.join(content)+"\n"
        return line

    def annotate_variants(self, allvariants):
        """
        Retrieve and merge information about the variants filtered

        allvariants (Variants) : Current variants processed

        return (list) : A list of lines in CSV format and the corresponding VCF content

        """
        sortedKeys = sorted(allvariants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        vcfcorrespondingline = 0
        printedLine = ""
        filteredvcfcontent = ""
        #Body
        for key in sortedKeys:
            endPosition = allvariants.variants[key]["position"] + len(allvariants.variants[key]["ref"]) - 1
            gffIds = self.annotations.map_GFFid_VCFpos(allvariants.variants[key]["chromosome"], 
            allvariants.variants[key]["position"], 
            endPosition)

            allvariants.variants[key]["features"] = gffIds
            self.define_categories(allvariants, key)
            if len(allvariants.variants[key]["features"]) == 0 and self.intergenicRegions is False:
                continue #No feature-surrounded variants
            allvariants.variants[key]["refwindow"] = self.fasta.window_sequence(allvariants.variants[key]["chromosome"],
            allvariants.variants[key]["position"], allvariants.variants[key]["ref"],
            Config.options["nuclWindowBefore"], Config.options["nuclWindowAfter"])

            self.load_features(allvariants, key)
            #Each alt has its line
            for altIndex, alt in enumerate(allvariants.variants[key]["alts"]):
                if allvariants.variants[key]["types"][altIndex] == "ambiguous":
                    continue #The VAF or RAF are unknown (low depth) or in between the minimum and maximum set
                printedLine += self.annotate_variant(allvariants, key, altIndex)
                self.print_log(allvariants, key)
                if allvariants.variants[key]["vcfline"] > vcfcorrespondingline:#Do not print the same line twice
                    vcfcorrespondingline = allvariants.variants[key]["vcfline"]
                    filteredvcfcontent += allvariants.vcf.vcfchunk[vcfcorrespondingline-1]
        return [printedLine, filteredvcfcontent]
    
    def get_variants_by_chunk(self, chunknumber):
        """
        Process variants for a particular chunk of the total set of variants

        chunknumber (int) : The rank of the chunk to be processed

        return (list) : The ordered sample names of the variants processed and their annotation

        """
        start_time = time.time()
        vcfchunk = Vcfchunk(chunknumber)
        vcfchunk.read_vcf()
        log = "Variants read in %s seconds"%(round(time.time()- start_time))

        start_time = time.time()
        allvariants = Variants(vcfchunk, Config.options["minSamplesDiff"])
        allvariants.process_variants(self.group1, self.group2)
        log += "\nVariants analysed in %s seconds"%(round(time.time()- start_time))

        start_time = time.time()
        variantAnnot = self.annotate_variants(allvariants)
        log += "\nVariant annotations added in %s seconds"%(round(time.time() - start_time))

        return [allvariants.samples, variantAnnot, log]


    def get_variants(self, outFileName):
        """
        Load, process and save all variants in a custom output file

        outFileName (str) : The name of the output file to save variants to

        """
        start_time = time.time()
        self.vcf = Vcf
        self.vcf.read_vcf_header(Config.options["vcf"])
        self.vcf.count_lines()
        self.chunks = max(math.ceil((self.vcf.totlines-self.vcf.headerlinenumber)/Config.options["chunksize"]), Config.options["ncores"])
        intervals = self.vcf.define_intervals(self.chunks)
        Config.verbose_print("Variant metadata read in %s seconds"%(round(time.time() - start_time)))
        log = ("Will use %s chunks "
                "each with %s variants "
                "for a total of %s processed variants.")%(self.chunks, 
                                                          "/".join((str(interval) for interval in intervals)), 
                                                          self.vcf.totlines-self.vcf.headerlinenumber)
        Config.verbose_print(log)
        
        csv = outFileName+".csv"
        filteredvcf = outFileName+".vcf" if self.outputVcf else None
        
        headerprinted = False
        max_chunk_set = min(self.chunks, Config.options["ncores"])
        chunk_set = 0
        while chunk_set < self.chunks:
            incr = max_chunk_set if chunk_set + max_chunk_set <= self.chunks else self.chunks - chunk_set
            chunk_set += incr
            
            Config.verbose_print("Processing chunks %s to %s..."%(chunk_set-incr+1, chunk_set))
            pool = Pool(Config.options["ncores"])
            results = pool.map(self.get_variants_by_chunk, range(chunk_set-incr+1,chunk_set+1))
            Config.verbose_print("\n".join(res[2] for res in results))
            self.samples = results[0][0]
            if not headerprinted:
                csvout, vcfout = self.print_header()
            else:
                csvout = ""
                vcfout = ""
            csvout += "".join(res[1][0] for res in results)
            vcfout += "".join(res[1][1] for res in results)
            
            if not headerprinted:
                opentype = 'w'
                headerprinted = True
            else:
                opentype = 'a'
            with open(csv, opentype) as f:
                f.write(csvout)
            if self.outputVcf:
                with open(filteredvcf, opentype) as f:
                    f.write(vcfout)
            
            del pool, results, csvout, vcfout
        
    def get_all_groups(self):
        """Define the groups to use to compare variants with each other"""
        compare_families = True if self.comparison == "families" else False
        compare_lineages = True if self.comparison == "lineages" else False
        compare_selfself = True if self.comparison == "self" else False

        if not compare_families and not compare_lineages and not compare_selfself:
            self.group1 = self.families.samples
            self.group2 = self.families.samples
            self.get_variants(self.outFile)
        
        if compare_families:
            families = list(self.families.families.keys())
            long_names = True if max([len(samplename) for samplename in families]) > 30 else False
            if long_names:
                print("Sample names exceed the limit of 30 characters, will use index instead for file names", file = stderr)
            assert len(families) > 1, "Not enough families to make a group comparison"

            for family_id in range(0,len(families)-1):
                family1 = families[family_id]
                self.group1 = self.families.families[family1]
                for other_family_id in range(family_id+1, len(families)):
                    family2 = families[other_family_id]
                    self.group2 = self.families.families[family2]
                    print("Group comparison of families %s (id : %i) and %s (id : %i)"%(family1, family_id, family2, other_family_id))
                    outFileInfo = self.outFile + "-"+family1+"_with_"+family2 if not long_names else self.outFile + "-"+str(family_id)+"_with_"+str(other_family_id)
                    self.get_variants(outFileInfo)
        
        if compare_lineages:
            parents = list(self.families.parents.keys())
            long_names = True if max([len(samplename) for samplename in parents]) > 30 else False
            if long_names:
                print("Sample names exceed the limit of 30 characters, will use index instead for file names", file = stderr)
            assert len(parents) >= 1, "Not enough parents-offsprings to make a group comparison"

            for index,mates in enumerate(parents):
                if mates == "NA":
                    continue
                matesnames = mates.split()
                plural = ["s", " and "] if len(matesnames) == 2 else ["", ""]
                addingmates = "-and-" if len(matesnames) == 2 else ""
                mate1 = matesnames[0]
                mate2 = matesnames[1] if len(matesnames) == 2 else ""

                self.group1 = matesnames
                self.group2 = self.families.parents[mates]

                print("Group comparison of parent%s %s%s%s (id : %i) with their offspring"%(plural[0], mate1, plural[1], mate2, index))
                outFileInfo = self.outFile + "-"+mate1+addingmates+mate2+"_with_offspring" if not long_names else self.outFile + "-"+str(index)+"_with_offspring"
                self.get_variants(outFileInfo)
        
        if compare_selfself:
            families = list(self.families.families.keys())
            long_names = True if max([len(samplename) for samplename in families]) > 30 else False
            if long_names:
                print("Sample names exceed the limit of 30 characters, will use index instead for file names", file = stderr)
            assert len(families) > 0, "There are no families in the PED file"

            for family_id in range(0,len(families)):
                family = families[family_id]
                self.group1 = self.families.families[family]
                self.group2 = self.group1
                print("Self-self comparison of family %s (id : %i)"%(family, family_id))
                outFileInfo = self.outFile + "-"+family if not long_names else self.outFile + "-"+str(family_id)
                self.get_variants(outFileInfo)
        
