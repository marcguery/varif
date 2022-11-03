from sys import stderr
from .config import Config
from .vcf import Vcf
from .variants import Variants
from .families import Families
from .annotations import Annotations
from .fasta import Fasta
from .version import __version__


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

        """
        Config.set_options()
        self.check_arguments(Config.options)
        self.vcf = Vcf()
        self.vcf.read_vcf(Config.options["vcf"])
        self.families = Families()
        if Config.options["ped"] is not None:
            self.families.read_ped(Config.options["ped"])
        self.annotations=Annotations()
        self.annotations.load_annotations_from_GFF(Config.options['gff'])
        self.fasta=Fasta()
        self.fasta.load_data_from_FASTA(Config.options['fasta'])
        self.get_all_groups(
            comparison=Config.options["comparison"],fixed=Config.options["fixed"], 
            allVariants=Config.options["allVariants"],allRegions=Config.options["allRegions"],
            outFile=Config.options["outFile"], outputVcf=Config.options["outputVcf"])
    
    def check_arguments(self, arguments):
        """
        Verify the integrity of the arguments passed from the command line

        arguments (dict) : Name of the argument (key) and its value

        """
        if Config.options["version"]:
            print(__version__)
            raise SystemExit
        for arg in ["vcf", "gff", "fasta", "outFile"]:
            if arguments[arg] is None:
                raise NameError("Argument '%s' is required"%arg)
        if arguments["comparison"] is not None:
            allowed_values=["families", "lineages", "both"]
            if arguments["comparison"] not in allowed_values:
                raise ValueError("Argument 'comparison' can only be one of '%s'"%("', '".join(allowed_values)))
            if arguments["ped"] is None:
                raise NameError("Argument 'ped' is required when comparing groups ('%s' comparison)"%(arguments["comparison"]))
        for arg in ["nuclWindowBefore", "nuclWindowAfter"]:
            if arguments[arg] < 0:
                raise ValueError("DNA window '%s' must be above 0, not %s"%(arg, arguments[arg])) 
        for arg in ["protWindowBefore", "protWindowAfter"]:
            if arguments[arg] < 0:
                raise ValueError("Protein window '%s' must be above 0, not %s"%(arg, arguments[arg]))
        if len(arguments["outFile"]) > 30:
            raise ValueError("Name of the outfile must not exceed 30 characters")
        print("Running Varif %s with options \n%s"%(__version__, "\n".join(" : ".join([opt, str(Config.options[opt])]) for opt in Config.options)))
    
    def merge_CDS(self, gffId, chromosome):
        """
        Create the spliced CDS sequence from a genomic CDS sequence

        gffId (str) : Identifier of the CDS
        chromosome (str) : Name of the chromosome

        return (list) : Sequence of the merged CDS and rtelative position of the initial CDS

        """
        parentalid = self.annotations.annotations[gffId]["parents"][0]
        allCDS = [cds for cds in self.annotations.genesID[self.annotations.annotations[gffId]['masterid']]["CDS"] if self.annotations.annotations[cds]["parents"][0] == parentalid]
        allCDSsorted = sorted(allCDS, key = lambda x : self.annotations.annotations[x]["start"])
        genesequence = ""
        for cds in allCDSsorted:
            start=self.annotations.annotations[cds]["start"]
            end=self.annotations.annotations[cds]["end"]
            if cds == gffId:
                startCDS = len(genesequence)+1
            genesequence+=self.fasta.data[chromosome][start-1:end]
        assert len(genesequence)%3 == 0
        return [genesequence, startCDS]
    
    def get_aa_from_mutation(self, chromosome, position, reference, mutation, gffId):
        """
        Get the change of aminoacid induced by a mutation

        chromosome (str) : Chromosome where the mutation occurs
        position (int) : Position of the mutation
        reference (str) : Upper case sequence before mutation
        mutation (str) : Upper case sequence after mutation
        gffId (str) : ID of the CDS as in GFF3 file

        return (list) : Aminoacids before and after mutation

        """
        assert self.annotations.annotations[gffId]['annotation']=='CDS'
        if len(self.annotations.annotations[gffId]["parents"]) > 1:
            print("Ambiguous origin of CDS %s: one of %s. Taking the first one (%s) by default"%(gffId, ", ".join(self.annotations.annotations[gffId]["parents"]), self.annotations.annotations[gffId]["parents"][0]), file = stderr)
        genesequence, startCDS = self.merge_CDS(gffId, chromosome)
        start=self.annotations.annotations[gffId]["start"] - (startCDS -1)
        strand=self.annotations.annotations[gffId]['strand']
        windowBefore=2+Config.options["protWindowBefore"]*3
        windowAfter=2+Config.options["protWindowAfter"]*3

        startIndex = max(0, position - start)
        endIndex=min(len(genesequence), startIndex+len(reference))
        
        genesequencemut, startIndexMut, endIndexMut = self.fasta.insert_mutation(chromosome, genesequence, strand, start, position, reference, mutation)

        oldCDS, aaPos, phase = self.fasta.window_CDS(genesequence, strand, startIndex, endIndex, windowBefore, windowAfter)
        newCDS, newaaPos, newphase = self.fasta.window_CDS(genesequencemut, strand, startIndexMut, endIndexMut, windowBefore, windowAfter)
        oldProt=self.fasta.translate_CDS(oldCDS, strand, phase)
        newProt=self.fasta.translate_CDS(newCDS, strand, newphase)
        
        aaChanges=[[aaPos,int(len(genesequence)/3)], oldCDS, oldProt, [newaaPos,int(len(genesequencemut)/3)], newCDS, newProt]
        return aaChanges
    
    def load_features(self, allvariants, variantId):
        """
        Store amino acids obtained from CDS translations
        
        allvariants (Variants) : Current variants processed
        variantId (str) : Unique identifier of the variant

        """
        for gffId in allvariants.variants[variantId]["features"]:
            if self.annotations.annotations[gffId]['annotation']!='CDS':
                continue #Not a CDS feature
            allvariants.variants[variantId]["aaPosAlts"][gffId]=[]
            allvariants.variants[variantId]["aaAlts"][gffId]=[]
            allvariants.variants[variantId]["cdsAlts"][gffId]=[]
            for alt in allvariants.variants[variantId]["alts"]:
                aaDiff=self.get_aa_from_mutation(allvariants.variants[variantId]["chromosome"], 
                allvariants.variants[variantId]["position"], 
                allvariants.variants[variantId]["ref"], alt, gffId)
                allvariants.variants[variantId]["aaPosAlts"][gffId].append(aaDiff[3])
                allvariants.variants[variantId]["cdsAlts"][gffId].append(aaDiff[4])
                allvariants.variants[variantId]["aaAlts"][gffId].append(aaDiff[5])
            allvariants.variants[variantId]["aaPosRef"][gffId]=aaDiff[0]
            allvariants.variants[variantId]["cdsRef"][gffId]=aaDiff[1]
            allvariants.variants[variantId]["aaRef"][gffId]=aaDiff[2]
    
    def define_category(self, allvariants, variantId):
        """
        Determine of the variant is a SNP or INDEL and add genomic location

        allvariants (Variants) : Current variants processed
        variantId (str) : Unique identifier of the variant

        """
        if all(len(allvariants.variants[variantId]["ref"]) == len(alt) for alt in allvariants.variants[variantId]["alts"]):
            maincategory = "SNP"
        else:
            maincategory = "INDEL"
        subcategories = sorted(set([self.annotations.annotations[gffId]['annotation'] for gffId in allvariants.variants[variantId]["features"]]))
        allvariants.variants[variantId]["category"] = ",".join([maincategory]+[subcategory for subcategory in subcategories])

    def print_line(self, allvariants, variantId=None, altIndex=0, sep=";"):
        """
        Create a line for the main output of varif separating each alt's variant

        allvariants (Variants) : Current variants processed
        variantId (str) : The unique identifier of the variant
        altIndex (int) : The index of the alt of this variant
        sep (str) : The field separator

        return (str) : The line of the header or of each alt

        """
        if variantId is None:
            baseHeader=[
                "Chromosome", "Position", "Type",
                "Ref", "Alt", "CDSref", "CDSalt", "AAref", "AAalt",
                "Annotation","Proportions"]
            line=sep.join(baseHeader+allvariants.samples)+"\n"
        else:
            windowLenSum=len(allvariants.variants[variantId]["refwindow"][0])+len(allvariants.variants[variantId]["refwindow"][1])
            leftwindow=allvariants.variants[variantId]["refwindow"][0]+"|" if windowLenSum > 0 else ""
            rightwindow="|"+allvariants.variants[variantId]["refwindow"][1] if windowLenSum > 0 else ""
            content=[
                allvariants.variants[variantId]["chromosome"], str(allvariants.variants[variantId]["position"]), 
                allvariants.variants[variantId]["category"],
                leftwindow+allvariants.variants[variantId]["ref"]+rightwindow, 
                allvariants.variants[variantId]["alts"][altIndex],
                "NA", "NA", "NA","NA", "NA",
                ":".join(f'{props:03}' for props in allvariants.variants[variantId]["props"][altIndex])]
            content+=["NA"]*len(allvariants.samples)

            aaposref=[]
            cdsref=[]
            aaref=[]
            aaposalts=[]
            cdsalts=[]
            aaalts=[]
            for gffId in allvariants.variants[variantId]["cdsRef"]:
                aaposref.append("/".join(str(pos) for pos in allvariants.variants[variantId]["aaPosRef"][gffId]))
                cdsref.append(allvariants.variants[variantId]["cdsRef"][gffId]+" ("+gffId+")")
                aaref.append(allvariants.variants[variantId]["aaRef"][gffId])
                aaposalts.append("/".join(str(pos) for pos in allvariants.variants[variantId]["aaPosAlts"][gffId][altIndex]))
                cdsalts.append(allvariants.variants[variantId]["cdsAlts"][gffId][altIndex])
                aaalts.append(allvariants.variants[variantId]["aaAlts"][gffId][altIndex])
            
            if len(allvariants.variants[variantId]["cdsRef"]) > 0: #There is a feature at least
                content[5]=":".join(cdsref) #the same cds ref in different features
                content[6]=":".join(cdsalts) #different cds alt
                content[7]="):".join(" (".join(str(cont) for cont in couple) for couple in zip(aaref, aaposref))+")" #the same aa ref in different features
                content[8]="):".join(" (".join(str(cont) for cont in couple) for couple in zip(aaalts, aaposalts))+")" #different aa alt
            annotations=set(self.annotations.annotations[geneId]['description']+" ("+geneId+")" for geneId in allvariants.variants[variantId]["features"] if "gene" in self.annotations.annotations[geneId]['annotation'])
            content[9]=":".join(annotations) if annotations!=set() else content[9] #potentially different annotations
            for i, sample in enumerate(allvariants.samples):#samples alt ratios
                content[11+i]=str(allvariants.variants[variantId]["ratios"][sample][altIndex+1])
            line=sep.join(content)+"\n"
        return line
        
    def print_log(self, allvariants, variantId):
        """
        Print the log stored in a variant if not empty

        allvariants (Variants) : Current variants processed
        variantId (str) : The unique identifier of the variant

        """
        log = ""
        mainlog = allvariants.variants[variantId]["log"][0]
        varlogs=allvariants.variants[variantId]["log"][1]

        if mainlog != "":
            log = "Variant %s: %s"%(variantId, mainlog)
        if "".join(varlogs) != "": 
            pass
        if log != "":
            print(log, file = stderr)

    def get_all_variants(self, allvariants, fixed, allVariants, allRegions, csv, filteredvcf):
        """
        Retrieve and merge information about the variants filtered

        allvariants (Variants) : Current variants processed
        fixed (bool) : Show also fixed variants if True
        allVariants (bool) : Show all variants if True (even fixed)
        allRegions (bool) : Show not only CDS if True
        csv (str) : File path of the CSV (semicolon separated)
        filteredvcf (str) : File path of the filtered VCF

        """
        fixed=True if allVariants is True else fixed
        sortedKeys=sorted(allvariants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        #Headers
        printedLine=self.print_line(allvariants)
        filteredvcfcontent="".join(allvariants.vcf.vcffile[0:allvariants.vcf.headerlinenumber])
        vcfcorrespondingline=allvariants.vcf.headerlinenumber-1
        #Body
        for key in sortedKeys:
            endPosition = allvariants.variants[key]["position"] + len(allvariants.variants[key]["ref"]) - 1
            gffIds=self.annotations.map_GFFid_VCFpos(allvariants.variants[key]["chromosome"], 
            allvariants.variants[key]["position"], 
            endPosition)

            allvariants.variants[key]["features"] = gffIds
            self.define_category(allvariants, key)
            if len(allvariants.variants[key]["features"]) == 0 and allRegions is False:
                continue #No feature-surrounded variants
            allvariants.variants[key]["refwindow"] = self.fasta.window_sequence(allvariants.variants[key]["chromosome"],
            allvariants.variants[key]["position"], allvariants.variants[key]["ref"],
            Config.options["nuclWindowBefore"], Config.options["nuclWindowAfter"])

            self.load_features(allvariants, key)
            #Each alt has its line
            for altIndex, alt in enumerate(allvariants.variants[key]["alts"]):
                if allvariants.variants[key]["types"][altIndex] == "ambiguous" and allVariants is False:
                    continue #The VAF or RAF are unknown (low depth) or in between the minimum and maximum set
                elif allvariants.variants[key]["types"][altIndex] == "fixed" and fixed is False:
                    continue #This variant is fixed (all are ref or all are alt)
                printedLine+=self.print_line(allvariants, key, altIndex)
                self.print_log(allvariants, key)
                if allvariants.variants[key]["vcfline"] > vcfcorrespondingline:#Do not print the same line twice
                    vcfcorrespondingline=allvariants.variants[key]["vcfline"]
                    filteredvcfcontent+=allvariants.vcf.vcffile[vcfcorrespondingline-1]
        with open(csv, 'w') as f:
            f.write(printedLine)
        if filteredvcf is not None:
            with open(filteredvcf, 'w') as f:
                f.write(filteredvcfcontent)
        
    def get_all_groups(self, comparison, fixed, allVariants, allRegions, outFile, outputVcf):
        """
        Get the processed variants for each combination of groups to be compared and output them

        comparison (str) : Group samples by family ('families'), lineage ('lineages'), both ('both') or none
        fixed (bool) : Show also fixed variants if True
        allVariants (bool) : Show all variants if True (even fixed)
        allRegions (bool) : Show not only CDS if True
        outFile (str) : Name (without extension) of the output file that will be appended with the group information
        outputVcf (bool) : Whether to ouput the corresponding VCF for each combination of groups compared

        """
        compare_families = True if comparison in ["both", "families"] else False
        compare_lineages = True if comparison in ["both", "lineages"] else False

        if not compare_families and not compare_lineages:
            allvariants=Variants(self.vcf)
            allvariants.process_variants()
            csv = outFile+".csv"
            filteredvcf = outFile+".vcf" if outputVcf else None
            
            self.get_all_variants(allvariants, fixed, allVariants, allRegions, csv, filteredvcf)
        
        if compare_families:
            families = list(self.families.families.keys())
            long_names = True if max([len(samplename) for samplename in families]) > 30 else False
            if long_names:
                print("Sample names exceed the limit of 30 characters, will use index instead for file names", file = stderr)
            assert len(families) > 1, "Not enough families to make a group comparison"

            for family_id in range(0,len(families)-1):
                family1 = families[family_id]
                group1 = self.families.families[family1]
                for other_family_id in range(family_id+1, len(families)):
                    family2 = families[other_family_id]
                    group2 = self.families.families[family2]

                    print("Group comparison of families %s (id : %i) and %s (id: %i)"%(family1, family_id, family2, other_family_id))
                    allvariants=Variants(self.vcf)
                    allvariants.process_variants(group1, group2)
                    outFileInfo = outFile + "-"+family1+"_with_"+family2 if not long_names else outFile + "-"+str(family_id)+"_with_"+str(other_family_id)
                    csv = outFileInfo+".csv"
                    filteredvcf = outFileInfo+".vcf" if outputVcf else None

                    self.get_all_variants(allvariants, fixed, allVariants, allRegions, csv, filteredvcf)
        
        if compare_lineages:
            parents = list(self.families.parents.keys())
            long_names = True if max([len(samplename) for samplename in parents]) > 30 else False
            if long_names:
                print("Sample names exceed the limit of 30 characters, will use index instead for file names", file = stderr)
            assert len(parents) >= 1, "Not enough parents-offsprings to make a group comparison"

            for index,mates in enumerate(parents):
                if mates == "NA":
                    continue
                matesnames=mates.split()
                plural = ["s", " and "] if len(matesnames) == 2 else ["", ""]
                addingmates = "-and-" if len(matesnames) == 2 else ""
                mate1 = matesnames[0]
                mate2 = matesnames[1] if len(matesnames) == 2 else ""

                group1 = matesnames
                group2 = self.families.parents[mates]

                print("Group comparison of parent%s %s%s%s (id : %i) with their offspring"%(plural[0], mate1, plural[1], mate2, index))
                allvariants=Variants(self.vcf)
                allvariants.process_variants(group1, group2)
                outFileInfo = outFile + "-"+mate1+addingmates+mate2+"_with_offspring" if not long_names else outFile + "-"+str(index)+"_with_offspring"
                csv = outFileInfo+".csv"
                filteredvcf = outFileInfo+".vcf" if outputVcf else None

                self.get_all_variants(allvariants, fixed, allVariants, allRegions, csv, filteredvcf)
