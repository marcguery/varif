import math
from sys import stderr
from .annotation import Annotation

class Annotations(object):
    """A bunch of annotations of a GFF3 file."""
    def __init__(self):
        """
        annotations (dict) : Key (ID of Annotation), value (Info of Annotation)
        genesID (dict) : Key (ID of gene annotation), value
                            (type of subfeatures: mRNA, exon, CDS...) containing their IDs
        positions (dict) : Key (Chromosome), value (Features coordinates and ID)
        index (dict) : For each chromosome (key), the value is a list with 
                            first element which is first coordinate of features,
                            second element which is a dict of features (value) with same last coordinate (key)
        ranks (list) : Ranks of fields order as in GFF3 specs

        """
        self.annotations = {}
        self.genesID = {}
        self.positions = {}
        self.index = {}
        self.ranks = range(0,9)
    
    def index_gff(self):
        """
        Store GFF features by their first coordinate containing groups of the same last coordinate

        return (dict) : For each chromosome (key), the value is a list with 
                            first element which is first coordinate of features,
                            second element which is a dict of features (value) with same last coordinate (key)
        
        """
        features = {}
        for chromosome in self.positions:
            self.positions[chromosome] = sorted(self.positions[chromosome], key = lambda x : x[0])
            currentInd = 0
            beeingAddedInd = [0]
            features[chromosome] = [(self.positions[chromosome][0][0], {self.positions[chromosome][0][1]:[self.positions[chromosome][0][2]]})]
            for feature in self.positions[chromosome][1:]:
                currentInd += 1
                if feature[0] == features[chromosome][-1][0]:
                    if feature[1] in features[chromosome][-1][1]:
                        features[chromosome][-1][1][feature[1]].append(feature[2])
                    else:
                        features[chromosome][-1][1][feature[1]] = [feature[2]]
                else:
                    features[chromosome].append((feature[0], {feature[1]:[feature[2]]}))
                    addedInd = beeingAddedInd
                    beeingAddedInd = [0]
                    indToCheck = addedInd+list(range(addedInd[-1]+1, currentInd))
                    for i in indToCheck:
                        if self.positions[chromosome][i][1] >= features[chromosome][-1][0]:
                            if self.positions[chromosome][i][1] in features[chromosome][-1][1]:
                                features[chromosome][-1][1][self.positions[chromosome][i][1]].append(self.positions[chromosome][i][2])
                            else:
                                features[chromosome][-1][1][self.positions[chromosome][i][1]] = [self.positions[chromosome][i][2]]
                            beeingAddedInd.append(i)
        features[chromosome] = sorted(features[chromosome], key = lambda x : x[0])
        return features
    
    def get_CDS_from_same_parent(self, gffId):
        """
        Retrieve all the CDS sequences with the same parent

        gffiId (str) : ID of the CDS

        return (list) : IDs of CDS with the same parent and their 0-based coordinate ranges
        
        """
        parentalid = self.annotations[gffId]["parents"][0]
        allCDS = [cds for cds in self.genesID[self.annotations[gffId]['masterid']]["CDS"] if self.annotations[cds]["parents"][0] == parentalid]
        allCDSsorted = sorted(allCDS, key = lambda x : self.annotations[x]["start"])
        cdsCoords = []
        for cds in allCDSsorted:
            start = self.annotations[cds]["start"]
            end = self.annotations[cds]["end"]
            cdsCoords.append([start-1, end])
        return [allCDSsorted, cdsCoords]

    def load_annotations_from_GFF(self, gff):
        """
        Store annotations and create a list of feature coordinates from GFF3

        gff (str) : Path of GFF3 file
        
        """
        gfffile = open(gff, 'r').readlines()
        n = 0
        #Header
        while n < len(gfffile) and gfffile[n].startswith('##'):
            n += 1
        while n < len(gfffile):
            if gfffile[n].startswith('#'):
                n += 1
                continue
            annotation = Annotation(gfffile[n], self.ranks)
            if annotation.id is None:
                #Discards entries without ID
                n += 1
                continue
            duplicatedidnumber = 1
            newannotation = annotation.id
            while newannotation in self.annotations and duplicatedidnumber < 100:
                print("There is a duplicate in the GFF file, ID: %s"%newannotation, file = stderr)
                print("Automatically assigning a new ID...", file = stderr)
                newannotation = annotation.id+".dupl."+str(duplicatedidnumber)
                duplicatedidnumber += 1
            annotation.id = newannotation
            #Filling descriptions with those of parents
            if annotation.parents != []:
                descs = []
                mID = None
                for parent in annotation.parents:
                    if parent in self.annotations:
                        descs.append(self.annotations[parent]['description'])
                        if mID is not None and self.annotations[parent]['masterid'] != mID:
                            print("Too much parents (%s, %s) for (%s)"%(mID[0], self.annotations[parent]['masterid'], annotation.id), file = stderr)
                        else:
                            mID = self.annotations[parent]['masterid']
                    else:
                        print("Orphan entry (%s) in the GFF file missing parent (%s)"%(annotation.id, parent), file = stderr)
                        raise(Exception)
                annotation.description = ",".join(descs)
                annotation.masterid = mID
                self.genesID[mID].setdefault(annotation.annotation, []).append(annotation.id)
            else:
                self.genesID[annotation.id] = {}
            #Storing annotations
            self.annotations[annotation.id] = {
                'chromosome':annotation.chromosome,
                'start':annotation.start,
                'end':annotation.end,
                'strand':annotation.strand,
                'phase':annotation.phase,
                'annotation':annotation.annotation,
                'parents':annotation.parents,
                'masterid':annotation.masterid,
                'description':annotation.description
            }
            #Generating feature coordinates
            if annotation.chromosome in self.positions:
                self.positions[annotation.chromosome].append([int(annotation.start), int(annotation.end), annotation.id])
            else:
                self.positions[annotation.chromosome] = [[int(annotation.start), int(annotation.end), annotation.id]]
            n += 1
        #Feature coordinates need to be sorted by start then end for variant mapping
        self.index = self.index_gff()
    
    def map_GFFid_VCFpos(self, chromosome, startPosition, endPosition):
        """
        Retrieve all GFF features surrounding the variant

        chromosome (str) : Name of chromosome
        startPosition (int) : First position of variant
        endPosition (int) : Last position of variant

        return (list) : Identifiers of GFF features retrieved
        
        """
        assert endPosition >= startPosition
        #This is a half-cut search
        finished = False
        #Maximum index of the feature to search
        maxValue = len(self.index[chromosome])-2
        #Minimum index of the feature to search
        minValue = 0
        #Half of those
        currindex = math.floor(maxValue+minValue/2)
        identifiers = []
        #Finding a random feature
        while not finished:
            #Minimum position of the feature checked
            mini = self.index[chromosome][currindex][0]
            #Maximum position of the feature checked
            maxi = max(list(self.index[chromosome][currindex][1].keys()))
            #The feature surrounds the variant
            if startPosition >= mini and startPosition <= maxi:
                finished = True
            #The variant is in upper half of the cut
            elif startPosition > maxi and currindex < maxValue:
                minValue = currindex+1
                currindex = math.ceil((maxValue+currindex)/2)
            #The variant is in lower half of the cut
            elif startPosition < mini and currindex > minValue:
                maxValue = currindex
                currindex = math.floor((currindex+minValue)/2)
            #The variant is not surrouded by a feature
            else:
                finished = True
        while currindex < len(self.index[chromosome])-1 and startPosition <= maxi and endPosition >= mini:
            for endFeature in self.index[chromosome][currindex][1]:
                if startPosition <= endFeature:
                    identifiers.extend(self.index[chromosome][currindex][1][endFeature])
            currindex += 1
            mini = self.index[chromosome][currindex][0]
            maxi = max(list(self.index[chromosome][currindex][1].keys()))
        return identifiers
