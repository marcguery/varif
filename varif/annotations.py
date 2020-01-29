from .annotation import Annotation

class Annotations(object):
    """A bunch of annotations of a GFF3 file."""
    def __init__(self):
        """
        annotations (dict) : Key (ID of Annotation), value (Info of Annotation)
        positions (dict) : Key (Chromosome), value (Features coordinates and ID)
        ranks (list) : Ranks of fileds order as in GFF3 specs

        """
        self.annotations={}
        self.positions={}
        self.index={}
        self.ranks=range(0,9)
    
    def index_gff(self):
        features={}
        for chromosome in self.positions:
            self.positions[chromosome]=sorted(self.positions[chromosome], key=lambda x : x[0])
            currentInd=0
            beeingAddedInd=[0]
            features[chromosome]=[(self.positions[chromosome][0][0], {self.positions[chromosome][0][1]:[self.positions[chromosome][0][2]]})]
            for feature in self.positions[chromosome][1:]:
                currentInd+=1
                if feature[0]==features[chromosome][-1][0]:
                    if feature[1] in features[chromosome][-1][1]:
                        features[chromosome][-1][1][feature[1]].append(feature[2])
                    else:
                        features[chromosome][-1][1][feature[1]]=[feature[2]]
                else:
                    features[chromosome].append((feature[0], {feature[1]:[feature[2]]}))
                    addedInd=beeingAddedInd
                    beeingAddedInd=[0]
                    indToCheck=addedInd+list(range(addedInd[-1]+1, currentInd))
                    for i in indToCheck:
                        if self.positions[chromosome][i][1] >= features[chromosome][-1][0]:
                            if self.positions[chromosome][i][1] in features[chromosome][-1][1]:
                                features[chromosome][-1][1][self.positions[chromosome][i][1]].append(self.positions[chromosome][i][2])
                            else:
                                features[chromosome][-1][1][self.positions[chromosome][i][1]]=[self.positions[chromosome][i][2]]
                            beeingAddedInd.append(i)
        features[chromosome]=sorted(features[chromosome], key=lambda x : x[0])
        return features

    def load_annotations_from_GFF(self, gff):
        """
        Store annotations and create a list of feature coordinates from GFF3

        gff (str) : Path of GFF3 file
        """
        gfffile=open(gff, 'r').readlines()
        n=0
        line=gfffile[n]
        #Header
        while line.startswith('##') and n < len(gfffile):
            line=gfffile[n]
            n+=1
        annotationWithoutDesc=[]
        while n < len(gfffile):            
            annotation=Annotation(gfffile[n], self.ranks)
            if annotation.id in self.annotations:
                print("There is a duplicate in the GFF file, ID : %s"%annotation.id)
                raise(Exception)
            #Filling descriptions with those of parents
            if annotation.description is None:
                descs=[]
                for parent in annotation.parents:
                    if parent in self.annotations:
                        descs.append(self.annotations[parent]['description'])
                    else:
                        descs=[]
                        print(annotation.parents)
                        annotationWithoutDesc.append(annotation.id)
                        break
                annotation.description=",".join(descs)
            #Storing annotations
            self.annotations[annotation.id]={
                'chromosome':annotation.chromosome,
                'start':annotation.start,
                'end':annotation.end,
                'strand':annotation.strand,
                'phase':annotation.phase,
                'annotation':annotation.annotation,
                'parents':annotation.parents,
                'description':annotation.description
            }
            #Generating feature coordinates
            if annotation.chromosome in self.positions:
                self.positions[annotation.chromosome].append([int(annotation.start), int(annotation.end), annotation.id])
            else:
                self.positions[annotation.chromosome]=[[int(annotation.start), int(annotation.end), annotation.id]]
            n+=1
        if len(annotationWithoutDesc) > 0:
            print("There were unresolved descriptions : \n%s"%", ".join(annotationWithoutDesc))
        #Feature coordinates need to be sorted by start then end for variant mapping
        self.index=self.index_gff()
