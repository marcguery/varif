from annotation import Annotation

class Annotations(object):
    def __init__(self):
        self.annotations={}
        self.positions={}
        self.ranks=range(0,9)

    def load_annotations_from_GFF(self, gff):
        gfffile=open(gff).readlines()
        n=0
        line=gfffile[n]
        while line.startswith('##') and n < len(gfffile):
            line=gfffile[n]
            n+=1
        annotationWithoutDesc=[]
        while n < len(gfffile):            
            annotation=Annotation(gfffile[n], self.ranks)
            if annotation.id in self.annotations:
                print("There is a duplicate in the GFF file, ID : %s"%annotation.id)
                raise(Exception)
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
            self.annotations[annotation.id]={
                'chromosome':annotation.chromosome,
                'start':annotation.start,
                'end':annotation.end,
                'annotation':annotation.annotation,
                'parents':annotation.parents,
                'description':annotation.description
            }
            if annotation.chromosome in self.positions:
                self.positions[annotation.chromosome].append([int(annotation.start), int(annotation.end), annotation.id])
            else:
                self.positions[annotation.chromosome]=[[int(annotation.start), int(annotation.end), annotation.id]]
            n+=1
        if len(annotationWithoutDesc) > 0:
            print("There were unresolved descriptions : \n%s"%", ".join(annotationWithoutDesc))
        for chromosome in self.positions:
            self.positions[chromosome]=sorted(self.positions[chromosome], key=lambda x : x[0])
