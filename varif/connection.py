import math
from .config import Config
from .variants import Variants
from .annotations import Annotations
from .fasta import Fasta


class Connection(object):
    """Filter and annotate variants different or fixed among samples"""

    def __init__(self, **options):
        """
        Arguments:
        vcf (str) : File path of VCF
        gff (str) : File path of GFF3
        fasta (str) : File path of Fasta
        All files must have the same chromosome names

        variants (Variants) : Variants object
        annotations (Annotations) : Annotations object
        fasta (Fasta) : Fasta object
        """
        Config.set_options(options)
        self.variants=Variants()
        self.variants.load_variants_from_VCF(Config.options['vcf'])
        self.annotations=Annotations()
        self.annotations.load_annotations_from_GFF(Config.options['gff'])
        self.fasta=Fasta()
        self.fasta.load_data_from_FASTA(Config.options['fasta'])
        self.get_scores(
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
        #This is a half-cut sort
        finished=False
        #Maximum index of the feature to search
        maxValue=len(self.annotations.positions[chromosome])-1
        #Minimum index of the feature to search
        minValue=0
        #Half of those
        index=math.floor(maxValue-minValue/2)
        identifiers=[]
        while not finished:
            #Minimum position of the feature checked
            mini=self.annotations.positions[chromosome][index][0]
            #Maximum position of the feature checked
            maxi=self.annotations.positions[chromosome][index][1]
            #The feature surrounds the variant
            if position >= mini and position <= maxi:
                #Going upper to find all features above
                n=0
                while position >= mini and index < maxValue:
                    if position <= maxi:
                        identifiers.append(self.annotations.positions[chromosome][index][2])
                    index+=1
                    n+=1
                    mini=self.annotations.positions[chromosome][index][0]
                    maxi=self.annotations.positions[chromosome][index][1]
                #Going lower to find all features below
                index=index-n-1
                mini=self.annotations.positions[chromosome][index][0]
                maxi=self.annotations.positions[chromosome][index][1]
                while position <= maxi and index > minValue:
                    if position >= mini:
                        identifiers.append(self.annotations.positions[chromosome][index][2])
                    index-=1
                    mini=self.annotations.positions[chromosome][index][0]
                    maxi=self.annotations.positions[chromosome][index][1]
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
                finished=True
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
        endIndex=position+len(reference)-1
        #The number of bases to skip
        # knowing that the window is +2 -2 bases after mutation
        #Translation starts at start of feature
        if strand=="+":
            phase=(phase+startIndex-2+1-startCDS)%3
        #Translation starts at end of feature
        elif strand=="-":
            phase=(phase+endIndex+2-endCDS)%3
        oldCDS=self.fasta.data[chromosome][startIndex-2:endIndex+2]
        oldProt=self.fasta.translate_CDS(oldCDS, strand, phase)

        newCDS=self.fasta.data[chromosome][startIndex-2:startIndex]+mutation+self.fasta.data[chromosome][endIndex:endIndex+2]
        newProt=self.fasta.translate_CDS(newCDS, strand, phase)
        aaChanges=[oldProt, newProt]
        return aaChanges

    def get_scores(self, fixed, allVariants, allRegions, show, csv):
        """
        Retrieve and merge information about the variants filtered

        fixed (bool) : Show also fixed variants if True
        allVariants (bool) : Show all variants if True (even fixed)
        allRegions (bool) : Show not only CDS if True
        show (bool) : Print in stdout
        csv (str) : File path of the CSV (semicolon separated)

        """
        #Filter by score to know what to show
        code=0 if fixed is True else 1
        code=-math.inf if allVariants is True else code
        sortedKeys=sorted(self.variants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        #Headers
        printedLine="Chromsome; Position; Type; Ref; Alt; AAref; AAalt; Annotation; Score;"
        printedLine+=";".join(self.variants.samples)+"\n"
        #Body
        for key in sortedKeys:
            #Filter
            if not any(score>=code for score in self.variants.variants[key]["scores"]):
                continue
            chrom=key.split(":")[0]
            pos=int(key.split(":")[1].split(".")[0])
            gffIds=self.map_GFFid_VCFpos(chrom, pos)
            #No feature-surrounded variants
            if len(gffIds)==0 and allRegions is False:
                continue
            for gffId in gffIds:
                #No CDS feature
                if self.annotations.annotations[gffId]['annotation']!='CDS':
                    continue
                self.variants.variants[key]["aaAlts"][gffId]=[]
                for alt in self.variants.variants[key]["alts"]:
                    aaDiff=self.get_aa_from_mutation(chrom, pos, self.variants.variants[key]["ref"], alt, gffId)
                    self.variants.variants[key]["aaAlts"][gffId].append(aaDiff[1])
                self.variants.variants[key]["aaRef"][gffId]=aaDiff[0]
            #Each alt has its line
            for altIndex, alt in enumerate(self.variants.variants[key]["alts"]):
                #Chromosome
                printedLine+=key.split(":")[0]+";"
                #Position
                printedLine+=key.split(":")[1]+";"
                #Type
                printedLine+=self.variants.variants[key]["category"]+";"
                #Ref base
                printedLine+=self.variants.variants[key]["ref"]+";"
                #Alt base
                printedLine+=alt+";"
                #feature-surrounded variant
                if len(gffIds) > 0:
                    #Ref AA
                    printedLine+=":".join([self.variants.variants[key]["aaRef"][gffId]+" ("+gffId+")" for gffId in self.variants.variants[key]["aaRef"]])+";"
                    #Alt AA per CDS feature
                    printedLine+=":".join([self.variants.variants[key]["aaAlts"][gffId][altIndex] for gffId in self.variants.variants[key]["aaAlts"]])+";"
                    #Unique annotations found
                    printedLine+=":".join(set(self.annotations.annotations[gffId]['description'] for gffId in gffIds))+";"
                #No feature-surrounded variant
                else:
                    printedLine+="NA;"
                    printedLine+="NA;"
                    printedLine+="NA;"
                #Scores
                printedLine+=str(self.variants.variants[key]["scores"][altIndex])+";"
                #Sample alt ratios
                for sampIndex,sample in enumerate(self.variants.variants[key]["ratios"]):
                    printedLine+=str(self.variants.variants[key]["ratios"][sample][altIndex+1])
                    printedLine+=";" if sampIndex < len(self.variants.variants[key]["ratios"])-1 else "\n"
        if show is True:
            print(printedLine)
        with open(csv, 'w') as f:
            f.write(printedLine)
