import math
from variants import Variants
from annotations import Annotations


class Varif(object):

    def __init__(self, vcf, gff):
        self.variants=Variants()
        self.variants.load_variants_from_VCF(vcf)
        self.annotations=Annotations()
        self.annotations.load_annotations_from_GFF(gff)
    
    def map_GFFid_VCFpos(self, chromosome, position):
        finished=False
        maxValue=len(self.annotations.positions[chromosome])-1
        minValue=0
        index=math.floor(maxValue-minValue/2)
        while not finished:
            mini=self.annotations.positions[chromosome][index][0]
            maxi=self.annotations.positions[chromosome][index][1]
            if position >= mini and position <= maxi:
                finished=True
                identifier=self.annotations.positions[chromosome][index][2]
                return identifier
            elif position > maxi and index < maxValue:
                minValue=index+1
                index=math.ceil((maxValue+index)/2)
            elif position < mini and index > minValue:
                maxValue=index-1
                index=math.floor((index+minValue)/2)
            else:
                finished=True
        return

    def get_scores(self, fixed=False, allVariants=False, allRegions=False, show=False, csv="Variants.csv"):
        code=0 if fixed is True else 1
        code=-math.inf if allVariants is True else code
        sortedKeys=sorted(self.variants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        printedLine="Chromsome, Position, Type, Ref, Alt, Score, Annotation,"
        printedLine+=",".join(self.variants.samples)+"\n"
        for key in sortedKeys:
            if not any(score>=code for score in self.variants.variants[key]["scores"]):
                continue
            gffId=self.map_GFFid_VCFpos(key.split(":")[0], int(key.split(":")[1].split(".")[0]))
            if gffId is None and allRegions is False:
                continue
            printedLine+=key.split(":")[0]+","
            printedLine+=key.split(":")[1]+","
            printedLine+=self.variants.variants[key]["category"]+","
            printedLine+=self.variants.variants[key]["ref"]+","
            printedLine+=";".join(self.variants.variants[key]["alts"])+","
            printedLine+=";".join([str(score) for score in self.variants.variants[key]["scores"]])+","
            if gffId is not None:
                printedLine+=self.annotations.annotations[gffId]['description']+" ("+gffId+"),"
            else:
                printedLine+="NA,"
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

    
