import math
from variants import Variants
from annotations import Annotations
from fasta import Fasta


class Varif(object):

    def __init__(self, vcf, gff, fasta):
        self.variants=Variants()
        self.variants.load_variants_from_VCF(vcf)
        self.annotations=Annotations()
        self.annotations.load_annotations_from_GFF(gff)
        self.fasta=Fasta()
        self.fasta.load_data_from_FASTA(fasta)
    
    def map_GFFid_VCFpos(self, chromosome, position):
        finished=False
        maxValue=len(self.annotations.positions[chromosome])-1
        minValue=0
        index=math.floor(maxValue-minValue/2)
        identifiers=[]
        while not finished:
            mini=self.annotations.positions[chromosome][index][0]
            maxi=self.annotations.positions[chromosome][index][1]
            if position >= mini and position <= maxi:
                n=0
                while position >= mini and position <= maxi:
                    identifiers.append(self.annotations.positions[chromosome][index][2])
                    index+=1
                    n+=1
                    mini=self.annotations.positions[chromosome][index][0]
                    maxi=self.annotations.positions[chromosome][index][1]
                index=index-n-1
                mini=self.annotations.positions[chromosome][index][0]
                maxi=self.annotations.positions[chromosome][index][1]
                while position >= mini and position <= maxi:
                    identifiers.append(self.annotations.positions[chromosome][index][2])
                    index-=1
                    mini=self.annotations.positions[chromosome][index][0]
                    maxi=self.annotations.positions[chromosome][index][1]
                finished=True
            if position > maxi and index < maxValue:
                minValue=index+1
                index=math.ceil((maxValue+index)/2)
            elif position < mini and index > minValue:
                maxValue=index-1
                index=math.floor((index+minValue)/2)
            else:
                finished=True
        return identifiers
    
    def get_aa_from_mutation(self, chromosome, position, reference, mutation, gffId):
        assert self.annotations.annotations[gffId]['annotation']=='CDS'
        phase=int(self.annotations.annotations[gffId]['phase'])
        strand=self.annotations.annotations[gffId]['strand']
        startCDS=self.annotations.annotations[gffId]['start']
        endCDS=self.annotations.annotations[gffId]['end']
        startIndex=position-1
        endIndex=position+len(reference)-1
        if strand=="+":
            print(position)
            phase=(phase+startIndex-2+1-startCDS)%3
        elif strand=="-":
            phase=(phase+endIndex+2-endCDS)%3
        oldCDS=self.fasta.data[chromosome][startIndex-2:endIndex+2]
        oldProt=self.fasta.translate_CDS(oldCDS, strand, phase)

        newCDS=self.fasta.data[chromosome][startIndex-2:startIndex]+mutation+self.fasta.data[chromosome][endIndex:endIndex+2]
        newProt=self.fasta.translate_CDS(newCDS, strand, phase)
        aaChanges=[oldProt, newProt]
        return aaChanges

    def get_scores(self, fixed=False, allVariants=False, allRegions=False, show=False, csv="Variants.csv"):
        code=0 if fixed is True else 1
        code=-math.inf if allVariants is True else code
        sortedKeys=sorted(self.variants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        printedLine="Chromsome, Position, Type, Ref, Alt, AAref, AAalts, Score, Annotation,"
        printedLine+=",".join(self.variants.samples)+"\n"
        for key in sortedKeys:
            if not any(score>=code for score in self.variants.variants[key]["scores"]):
                continue
            chrom=key.split(":")[0]
            pos=int(key.split(":")[1].split(".")[0])
            gffIds=self.map_GFFid_VCFpos(chrom, pos)
            if len(gffIds)==0 and allRegions is False:
                continue
            for gffId in gffIds:
                if self.annotations.annotations[gffId]['annotation']!='CDS':
                    continue
                
                self.variants.variants[key]["aaAlts"][gffId]=[]

                for alt in self.variants.variants[key]["alts"]:
                    aaDiff=self.get_aa_from_mutation(chrom, pos, self.variants.variants[key]["ref"], alt, gffId)
                    self.variants.variants[key]["aaAlts"][gffId].append(aaDiff[1])
                self.variants.variants[key]["aaRef"][gffId]=aaDiff[0]
            printedLine+=key.split(":")[0]+","
            printedLine+=key.split(":")[1]+","
            printedLine+=self.variants.variants[key]["category"]+","
            printedLine+=self.variants.variants[key]["ref"]+","
            printedLine+=";".join(self.variants.variants[key]["alts"])+","
            if len(gffIds) > 0:
                printedLine+=":".join([self.variants.variants[key]["aaRef"][gffId]+" ("+gffId+")" for gffId in self.variants.variants[key]["aaRef"]])+","
                printedLine+=":".join([";".join(self.variants.variants[key]["aaAlts"][gffId])+" ("+gffId+")" for gffId in self.variants.variants[key]["aaAlts"]])+","
                printedLine+=":".join([self.annotations.annotations[gffId]['description']+" ("+gffId+")" for gffId in gffIds])+","
            else:
                printedLine+="NA,"
                printedLine+="NA,"
                printedLine+="NA,"
            printedLine+=";".join([str(score) for score in self.variants.variants[key]["scores"]])+","
            for index,sample in enumerate(self.variants.variants[key]["ratios"]):
                printedLine+=";".join([str(self.variants.variants[key]["ratios"][sample][rank]) for rank in range(1,len(self.variants.variants[key]["ratios"][sample]))])
                printedLine+="," if index < len(self.variants.variants[key]["ratios"])-1 else "\n"
        if show is True:
            print(printedLine)
        with open(csv, 'w') as f:
            f.write(printedLine)

        
    def get_id_from_position(self, position, gffreader, chromosome):
        """
        Retrieve information about an annotation level in which there is a variant

        position (int) : the position of the variant
        gffreader (read_gff) : the gff reader object
        chromosome (str) : the chromosome name in gff reader

        return (list) : ID, function, start/end positions or None if outside of a variant
        """
        finished=False
        maxValue=len(gffreader[chromosome]['positions'])-1
        minValue=0
        index=math.floor(maxValue-minValue/2)
        while not finished:
            mini=gffreader[chromosome]['positions'][index][1][0]
            maxi=gffreader[chromosome]['positions'][index][1][1]
            if position >= mini and position <= maxi:
                finished=True
                rank=0
                while gffreader[chromosome]['ids'][rank][0] != gffreader[chromosome]['positions'][index][0]:
                    rank+=1
                identifier=gffreader[chromosome]['ids'][rank][1]
                description=gffreader[chromosome]['descriptions'][rank][1]
                return [identifier, description, gffreader[chromosome]['positions'][index]]
            elif position > maxi and index < maxValue:
                minValue=index+1
                index=math.ceil((maxValue+index)/2)
            elif position < mini and index > minValue:
                maxValue=index-1
                index=math.floor((index+minValue)/2)
            else:
                finished=True

        return ["No ID", "No function", "No Positions"]

    
