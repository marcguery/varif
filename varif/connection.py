import math
import numpy as np
from sys import stderr
from .config import Config
from .variants import Variants
from .annotations import Annotations
from .fasta import Fasta


class Connection(object):
    """Filter and annotate variants different or fixed among samples"""

    def __init__(self):
        """
        variants (Variants) : Variants object
        annotations (Annotations) : Annotations object
        fasta (Fasta) : Fasta object
        """
        Config.set_options()
        self.check_arguments(Config.options)
        self.variants=Variants()
        self.variants.load_variants_from_VCF(Config.options['vcf'])
        self.annotations=Annotations()
        self.annotations.load_annotations_from_GFF(Config.options['gff'])
        self.fasta=Fasta()
        self.fasta.load_data_from_FASTA(Config.options['fasta'])
        self.get_all_alts(
            fixed=Config.options["fixed"], allVariants=Config.options["allVariants"],
            allRegions=Config.options["allRegions"],
            csv=Config.options["csv"], filteredvcf=Config.options["filteredvcf"])
    
    @property
    def version(self):
        return "0.2.0"
    
    def check_arguments(self, arguments):
        if Config.options["version"]:
            print(self.version)
            raise SystemExit
        for arg in ["vcf", "gff", "fasta"]:
            if arguments[arg] is None:
                raise NameError("Argument '%s' is required"%arg)
        for arg in ["nuclWindowBefore", "nuclWindowAfter"]:
            if arguments[arg] < 0 or arguments[arg] > 1000:
                raise ValueError("DNA window '%s' must be between 0 and 1000, not %s"%(arg, arguments[arg])) 
        for arg in ["protWindowBefore", "protWindowAfter"]:
            if arguments[arg] < 0 or arguments[arg] > 400:
                raise ValueError("Protein window '%s' must be between 0 and 400, not %s"%(arg, arguments[arg]))
        print("Running Varif %s with options \n%s"%(self.version, "\n".join(" : ".join([opt, str(Config.options[opt])]) for opt in Config.options)))
    
    def map_GFFid_VCFpos(self, chromosome, position):
        """
        Retrieve all GFF feature surrounding the variant

        chromosome (str) : Name of chromosome
        position (int) : First position of variant

        return (list) : Identifiers of GFF features retrieved
        """
        #This is a half-cut search
        finished=False
        #Maximum index of the feature to search
        maxValue=len(self.annotations.index[chromosome])-2
        #Minimum index of the feature to search
        minValue=0
        #Half of those
        index=math.floor(maxValue+minValue/2)
        identifiers=[]
        #Finding a random feature
        while not finished:
            #Minimum position of the feature checked
            mini=self.annotations.index[chromosome][index][0]
            #Maximum position of the feature checked
            maxi=self.annotations.index[chromosome][index+1][0]
            #The feature surrounds the variant
            if position >= mini and position <= maxi:
                finished=True
            #The variant is in upper half of the cut
            elif position > maxi and index < maxValue:
                minValue=index+1
                index=math.ceil((maxValue+index)/2)
            #The variant is in lower half of the cut
            elif position < mini and index > minValue:
                maxValue=index
                index=math.floor((index+minValue)/2)
            #The variant is not surrouded by a feature
            else:
                return identifiers
        identPerPos=[self.annotations.index[chromosome][index][1][key] for key in self.annotations.index[chromosome][index][1] if key>=position]
        identifiers=[identifier for position in identPerPos for identifier in position]
        return identifiers
    
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
        phase=int(self.annotations.annotations[gffId]['phase'])
        strand=self.annotations.annotations[gffId]['strand']
        startCDS=self.annotations.annotations[gffId]['start']
        endCDS=self.annotations.annotations[gffId]['end']
        startIndex=position-1
        endIndex=startIndex+len(reference)
        windowBefore=2+Config.options["protWindowBefore"]*3
        windowAfter=2+Config.options["protWindowAfter"]*3

        #The number of bases to skip
        # knowing that the window is +/- windows bases after mutation
        #Translation starts at start of feature
        if strand=="+":
            shiftBefore=min(startIndex-startCDS+1, windowBefore)
            shiftAfter=min(endCDS-endIndex, windowAfter)
            phase=-(startIndex-shiftBefore-(startCDS-1+phase))%3
        #Translation starts at end of feature
        elif strand=="-":
            shiftBefore=min(startIndex-startCDS+1, windowAfter)
            shiftAfter=min(endCDS-endIndex, windowBefore)
            phase=(endIndex+shiftAfter-phase-endCDS)%3
        
        oldCDS=self.fasta.data[chromosome][startIndex-shiftBefore:endIndex+shiftAfter]
        oldProt=self.fasta.translate_CDS(oldCDS, strand, phase)

        newCDS=self.fasta.data[chromosome][startIndex-shiftBefore:startIndex]+mutation+self.fasta.data[chromosome][endIndex:endIndex+shiftAfter]
        newProt=self.fasta.translate_CDS(newCDS, strand, phase)
        aaChanges=[oldCDS, oldProt, newCDS, newProt]
        return aaChanges

    def print_line(self, variantId=None, altIndex=0, sep=";"):
        """
        Create a line for the final output of varif using each alt's variant

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
            line=sep.join(baseHeader+self.variants.samples)+"\n"
        else:
            windowLenSum=len(self.variants.variants[variantId]["refwindow"][0])+len(self.variants.variants[variantId]["refwindow"][1])
            leftwindow=self.variants.variants[variantId]["refwindow"][0]+"|" if windowLenSum > 0 else ""
            rightwindow="|"+self.variants.variants[variantId]["refwindow"][1] if windowLenSum > 0 else ""
            content=[
                self.variants.variants[variantId]["chromosome"], str(self.variants.variants[variantId]["position"]), 
                self.variants.variants[variantId]["category"],
                leftwindow+self.variants.variants[variantId]["ref"]+rightwindow, 
                self.variants.variants[variantId]["alts"][altIndex],
                "NA", "NA", "NA","NA", "NA",
                ":".join(f'{props:03}' for props in self.variants.variants[variantId]["props"][altIndex])]
            content+=["NA"]*len(self.variants.samples)

            cdsref=[]
            aaref=[]
            cdsalts=[]
            aaalts=[]
            for gffId in self.variants.variants[variantId]["cdsRef"]:
                cdsref.append(self.variants.variants[variantId]["cdsRef"][gffId]+" ("+gffId+")")
                aaref.append(self.variants.variants[variantId]["aaRef"][gffId])
                cdsalts.append(self.variants.variants[variantId]["cdsAlts"][gffId][altIndex])
                aaalts.append(self.variants.variants[variantId]["aaAlts"][gffId][altIndex])
            
            if len(self.variants.variants[variantId]["cdsRef"]) > 0: #There is a feature at least
                content[5]=":".join(cdsref) #the same cds ref in different features
                content[6]=":".join(cdsalts) #different cds alt
                content[7]=":".join(aaref) #the same aa ref in different features
                content[8]=":".join(aaalts) #different aa alt
            annotations=set(self.annotations.annotations[geneId]['description']+" ("+geneId+")" for geneId in self.variants.variants[variantId]["features"] if "gene" in self.annotations.annotations[geneId]['annotation'])
            content[9]=":".join(annotations) if annotations!=set() else content[9] #potentially different annotations
            for i, sample in enumerate(self.variants.samples):#samples alt ratios
                content[11+i]=str(self.variants.variants[variantId]["ratios"][sample][altIndex+1])
            line=sep.join(content)+"\n"
        return line
    
    def window_sequence(self, variantId):
        """
        Extract DNA sequence before and after the variant

        variantId (str) : Unique identifier of the variant

        """
        chr = self.variants.variants[variantId]["chromosome"]
        pos = self.variants.variants[variantId]["position"]
        seqLength = len(self.variants.variants[variantId]["ref"])
        posBefore = pos - Config.options["nuclWindowBefore"]
        posAfter = pos + Config.options["nuclWindowAfter"]
        if posBefore < 1 or posAfter + seqLength > len(self.fasta.data[chr]):
            print("Window too wide for variant %s (%s:%s)"%(variantId, chr, pos), file=stderr)
            posBefore = 1 if posBefore < 1 else posBefore
            posAfter = -2 - seqLength if posAfter + seqLength > len(self.fasta.data[chr]) else posAfter

        self.variants.variants[variantId]["refwindow"]=[self.fasta.data[chr][posBefore-1:pos-1],self.fasta.data[chr][pos+seqLength-1:posAfter+seqLength-1]]
        
    def load_features(self, variantId):
        """
        Store amino acids obtained from CDS translations

        variantId (str) : Unique identifier of the variant

        """
        for gffId in self.variants.variants[variantId]["features"]:
            if self.annotations.annotations[gffId]['annotation']!='CDS':
                continue #Not a CDS feature
            self.variants.variants[variantId]["aaAlts"][gffId]=[]
            self.variants.variants[variantId]["cdsAlts"][gffId]=[]
            for alt in self.variants.variants[variantId]["alts"]:
                aaDiff=self.get_aa_from_mutation(self.variants.variants[variantId]["chromosome"], 
                self.variants.variants[variantId]["position"], 
                self.variants.variants[variantId]["ref"], alt, gffId)
                self.variants.variants[variantId]["cdsAlts"][gffId].append(aaDiff[2])
                self.variants.variants[variantId]["aaAlts"][gffId].append(aaDiff[3])
            self.variants.variants[variantId]["aaRef"][gffId]=aaDiff[1]
            self.variants.variants[variantId]["cdsRef"][gffId]=aaDiff[0]

    def get_all_alts(self, fixed, allVariants, allRegions, csv, filteredvcf):
        """
        Retrieve and merge information about the variants filtered

        fixed (bool) : Show also fixed variants if True
        allVariants (bool) : Show all variants if True (even fixed)
        allRegions (bool) : Show not only CDS if True
        csv (str) : File path of the CSV (semicolon separated)
        filteredvcf (str) : File path of the filtered VCF


        """
        fixed=True if allVariants is True else fixed
        sortedKeys=sorted(self.variants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        #Headers
        printedLine=self.print_line()
        filteredvcfcontent="".join(self.variants.vcffile[0:self.variants.headerlinenumber])
        vcfcorrespondingline=self.variants.headerlinenumber-1
        #Body
        for key in sortedKeys:
            gffIds=self.map_GFFid_VCFpos(self.variants.variants[key]["chromosome"], self.variants.variants[key]["position"])
            self.variants.variants[key]["features"]=gffIds
            if len(self.variants.variants[key]["features"]) == 0 and allRegions is False:
                continue #No feature-surrounded variants
            self.window_sequence(key)
            self.load_features(key)
            #Each alt has its line
            for altIndex, alt in enumerate(self.variants.variants[key]["alts"]):
                if self.variants.variants[key]["types"][altIndex] == "ambiguous" and allVariants is False:
                    continue #The VAF or RAF are unknown (low depth) or in between the minimum and maximum set
                elif self.variants.variants[key]["types"][altIndex] == "fixed" and fixed is False:
                    continue #This variant is fixed (all are ref or all are alt)
                printedLine+=self.print_line(key, altIndex)
                if self.variants.variants[key]["vcfline"] > vcfcorrespondingline:#Do not print the same line twice
                    vcfcorrespondingline=self.variants.variants[key]["vcfline"]
                    filteredvcfcontent+=self.variants.vcffile[vcfcorrespondingline-1]
        if csv is None and filteredvcf is None:
            print("-----Filtered CSV file:-----")
            print(printedLine.strip("\n"))
            print("-----Filtered VCF file:-----")
            print(filteredvcfcontent.strip("\n"))
        if csv is not None:
            with open(csv, 'w') as f:
                f.write(printedLine)
        if filteredvcf is not None:
            with open(filteredvcf, 'w') as f:
                f.write(filteredvcfcontent)
