import math
from variants import Variants

class Varif(object):

    def __init__(self, vcf):
        self.variants=Variants()
        self.variants.load_variants_from_VCF(vcf)
        self.annotations={}

    def get_scores(self, csv="Variants.csv", show=False):
        sortedKeys=sorted(self.variants.variants, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        printedLine="Chromsome, Position, Type, Ref, Alt,"
        printedLine+=",".join(self.variants.samples)+","
        printedLine+="Score\n"
        for key in sortedKeys:
            if not any(score==0 for score in self.variants.variants[key]["scores"]):
                continue
            printedLine+=key.split(":")[0]+","
            printedLine+=key.split(":")[1]+","
            printedLine+=self.variants.variants[key]["category"]+","
            printedLine+=self.variants.variants[key]["ref"]+","
            printedLine+=";".join(self.variants.variants[key]["alts"])+","
            for sample in self.variants.variants[key]["ratios"]:
                printedLine+=";".join([str(self.variants.variants[key]["ratios"][sample][rank]) for rank in range(1,len(self.variants.variants[key]["ratios"][sample]))])
                printedLine+=","
            printedLine+=";".join([str(score) for score in self.variants.variants[key]["scores"]])+"\n"
        with open(csv, 'w') as f:
            f.write(printedLine)


    def read_gff(self, gff, annotation):
        """
        Reads a GFF file

        gff (str) : the gff path location
        annotation (str) : annotation to be taken to subset positions

        return (dict) : The positions (start and end), ids and functions
        of extracted information per chromosome
        """
        gffInfo={}
        with open(gff) as f:
            rank=0
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                lineSplitted=line.split("\t", 8)
                if lineSplitted[2]==annotation:
                    rank+=1
                    identifier=lineSplitted[-1].split(";")[0][3:]
                    description=lineSplitted[-1].split(";")[-1][12:].rstrip("\n") if lineSplitted[-1].split(";")[-1].startswith('description') else 'Unknown function'
                    rangePosition=[int(lineSplitted[3]), int(lineSplitted[4])+1]
                    if not lineSplitted[0] in gffInfo:
                        gffInfo[lineSplitted[0]]={}
                        gffInfo[lineSplitted[0]]['positions']=[(rank, rangePosition)]
                        gffInfo[lineSplitted[0]]['ids']=[(rank, identifier)]
                        gffInfo[lineSplitted[0]]['descriptions']=[(rank, description)]
                    else:
                        gffInfo[lineSplitted[0]]['positions'].append((rank, rangePosition))
                        gffInfo[lineSplitted[0]]['ids'].append((rank, identifier))
                        gffInfo[lineSplitted[0]]['descriptions'].append((rank, description))
        for chromosome in gffInfo:
            gffInfo[chromosome]['positions']=sorted(gffInfo[chromosome]['positions'], key=lambda x : x[1][0])
        return gffInfo

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

    
