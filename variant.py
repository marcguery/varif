import math
import numpy as np

class Variant(object):
    def __init__(self, vcfLine, ranks, samples, samplesRanks):
        vcfLine=vcfLine.split("\t")
        self.chromosome=vcfLine[ranks[0]]
        self.position=vcfLine[ranks[1]]
            
        self.ref=vcfLine[ranks[3]]
        self.alts=vcfLine[ranks[4]].split(",")
        formatSplitted=vcfLine[ranks[8]].split(":")
        for i in range(len(formatSplitted)):
            if formatSplitted[i] == "AD":
                adRank=i
                break
        counts={samples[i]:[int(ad) for ad in vcfLine[n].split(":")[adRank].strip("\n").split(",")] for i,n in enumerate(samplesRanks)}
        self.counts=counts
        self.ratios={}
        self.scores=[]
    
    @property
    def category(self):
        if all(len(self.ref)==len(alt) for alt in self.alts):
            return "SNP"
        else:
            return "INDEL"
    
    def calculate_ratios(self, mindepth=5):
        #Handling division by zero, when there is no ref
        for sample in self.counts:
            if sum(self.counts[sample]) < mindepth:
                self.ratios[sample]=[math.nan]*len(self.counts[sample])
            else:
                self.ratios[sample]=[ad/sum(self.counts[sample]) for ad in self.counts[sample]]
        return self.ratios

    def scores_from_ratios(self, maxprop=0.1, minprop=None):
        assert self.ratios != {}
        minprop=1-maxprop if minprop is None else minprop

        assert 0 < maxprop < 0.5 and 0.5 < minprop < 1
        for rank in range(1,len(self.alts)+1):
            ratioRank=[self.ratios[sample][rank] for sample in self.ratios]
            #No sample is above depth
            if math.isnan(sum(ratioRank)):
                self.scores.append(math.nan)
                continue
            minratio=np.nanmin(ratioRank)
            maxratio=np.nanmax(ratioRank)
            #fixedVariants
            if minratio > minprop:
                self.scores.append(0)
            #Ambiguous variants
            elif minratio >= maxprop or maxratio <= minprop:
                self.scores.append(-maxratio/minratio if minratio != 0 else -math.inf)
            #True Variants
            elif minratio < maxprop and maxratio > minprop:
                self.scores.append(maxratio/minratio if minratio != 0 else math.inf)
        return self.scores