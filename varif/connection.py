import math
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
        self.variants=Variants()
        self.variants.load_variants_from_VCF(Config.options['vcf'])
        self.annotations=Annotations()
        self.annotations.load_annotations_from_GFF(Config.options['gff'])
        self.fasta=Fasta()
        self.fasta.load_data_from_FASTA(Config.options['fasta'])
        self.get_all_alts(
            fixed=Config.options["fixed"], allVariants=Config.options["allVariants"],
            allRegions=Config.options["allRegions"], show=Config.options["show"],
            csv=Config.options["csv"])
    
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
        #The number of bases to skip
        # knowing that the window is +2 -2 bases after mutation
        #Translation starts at start of feature
        if strand=="+":
            phase=-(startIndex-2-(startCDS-1+phase))%3
        #Translation starts at end of feature
        elif strand=="-":
            phase=(endIndex+2-phase-endCDS)%3
        oldCDS=self.fasta.data[chromosome][startIndex-2:endIndex+2]
        oldProt=self.fasta.translate_CDS(oldCDS, strand, phase)

        newCDS=self.fasta.data[chromosome][startIndex-2:startIndex]+mutation+self.fasta.data[chromosome][endIndex:endIndex+2]
        newProt=self.fasta.translate_CDS(newCDS, strand, phase)
        aaChanges=[oldProt, newProt]
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
                "Ref", "Alt", "AAref", "AAalt",
                "Annotation","Score"]
            line=sep.join(baseHeader+self.variants.samples+["Con", "Mut", "Mix", "Und"])+"\n"
        else:
            content=[
                self.variants.variants[variantId]["chromosome"], str(self.variants.variants[variantId]["position"]), 
                self.variants.variants[variantId]["category"],
                self.variants.variants[variantId]["ref"], self.variants.variants[variantId]["alts"][altIndex],
                "NA", "NA", "NA",
                str(self.variants.variants[variantId]["scores"][altIndex])]
            content+=["NA"]*len(self.variants.samples)
            content+=["NA", "NA", "NA", "NA"]

            aaref=[]
            aaalts=[]
            for gffId in self.variants.variants[variantId]["aaRef"]:
                aaref.append(self.variants.variants[variantId]["aaRef"][gffId]+" ("+gffId+")")
                aaalts.append(self.variants.variants[variantId]["aaAlts"][gffId][altIndex])
            
            if len(self.variants.variants[variantId]["aaRef"]) > 0: #There is a feature at least
                content[5]=":".join(aaref) #the same aa ref in different features
                content[6]=":".join(aaalts) #different aa alt
            annotations=set(self.annotations.annotations[geneId]['description']+" ("+geneId+")" for geneId in self.variants.variants[variantId]["features"] if self.annotations.annotations[geneId]['annotation']=="gene")
            content[7]=":".join(annotations) if annotations!=set() else content[7] #potentially different annotations
            groups={0:[], 1:[], 2:[], 3:[]}
            for i, sample in enumerate(self.variants.samples):#samples alt ratios
                content[9+i]=str(self.variants.variants[variantId]["ratios"][sample][altIndex+1])
                groupNumber=self.variants.variants[variantId]["groups"][sample][altIndex+1]
                groups[groupNumber].append(sample)
            for i in groups:
                content[9+len(self.variants.samples)+i]=":".join(groups[i])
            line=sep.join(content)+"\n"
        return line

    def load_features(self, variantId, gffId):
        """
        Store amino acids obtained from CDS translations of a CDS feature

        variantId (str) : Unique identifier of the variant
        gffId (str) : Unique identifier of the CDS feature to translate

        """
        self.variants.variants[variantId]["aaAlts"][gffId]=[]
        for alt in self.variants.variants[variantId]["alts"]:
            aaDiff=self.get_aa_from_mutation(self.variants.variants[variantId]["chromosome"], 
            self.variants.variants[variantId]["position"], 
            self.variants.variants[variantId]["ref"], alt, gffId)
            self.variants.variants[variantId]["aaAlts"][gffId].append(aaDiff[1])
        self.variants.variants[variantId]["aaRef"][gffId]=aaDiff[0]

    def get_all_alts(self, fixed, allVariants, allRegions, show, csv):
        """
        Retrieve and merge information about the variants filtered

        fixed (bool) : Show also fixed variants if True
        allVariants (bool) : Show all variants if True (even fixed)
        allRegions (bool) : Show not only CDS if True
        show (bool) : Print in stdout
        csv (str) : File path of the CSV (semicolon separated)

        """
        #Filter by score to know what to show
        code=-1 if fixed is True else 1
        code=-math.inf if allVariants is True else code
        sortedKeys=sorted(self.variants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        #Headers
        printedLine=self.print_line()
        #Body
        for key in sortedKeys:
            if not any(score > code for score in self.variants.variants[key]["scores"]):
                continue #Do not do anything if everyone is garbage
            gffIds=self.map_GFFid_VCFpos(self.variants.variants[key]["chromosome"], self.variants.variants[key]["position"])
            self.variants.variants[key]["features"]=gffIds
            if len(self.variants.variants[key]["features"])==0 and allRegions is False:
                continue #No feature-surrounded variants
            for gffId in self.variants.variants[key]["features"]:
                if self.annotations.annotations[gffId]['annotation']!='CDS':
                    continue #Not a CDS feature
                self.load_features(key, gffId)
            #Each alt has its line
            for altIndex, alt in enumerate(self.variants.variants[key]["alts"]):
                if self.variants.variants[key]["scores"][altIndex] <= code:
                    continue #The score is not good enough
                printedLine+=self.print_line(key, altIndex)
        if show is True:
            print(printedLine)
        with open(csv, 'w') as f:
            f.write(printedLine)
