import math
import numpy as np
from config import Config

class Variant(object):
    """A parsed line of a VCF, requiring AD for each sample."""
    def __init__(self, vcfLine, ranks, samples, samplesRanks):
        """
        Arguments:
        vcfLine (str) : Raw VCF line
        ranks (list) : Order of fields as in VCF specs
        sample (list) : Name of samples 
        samplesRanks (list) : Order of samples given in argument

        chromosome (str) : Name of chromosome
        position (str) : Position of variant
        ref (str) : Reference sequence of variant starting at position
        alts (list) : Alternate sequences at position
        counts (dict) : Key (sample name), value (AD count for ref and alts)
        ratios (dict) : Key (sample name), value (ratio of AD for ref and alts)
        scores (list) : Score of each alt
        maximum (int) : Score attributed to x/0
        minimum (int) : Score attributed to -x/0
        category (str) : Type of variant : 'SNP' or 'INDEL'

        """
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
        self.counts={samples[i]:[int(ad) for ad in vcfLine[n].split(":")[adRank].strip("\n").split(",")] for i,n in enumerate(samplesRanks)}
        self.ratios={}
        self.scores=[]
        self.maximum=Config.options["maximum"]
        self.minimum=Config.options["minimum"]
    
    @property
    def category(self):
        if all(len(self.ref)==len(alt) for alt in self.alts):
            return "SNP"
        else:
            return "INDEL"
    
    def calculate_ratios(self, mindepth):
        """
        Get the ratio for each alternate count of the variant

        mindepth (int) : Minmal depth for a sample to consider its counts
        To show every variant later, mindepth must be 0

        return (dict) : Ratios for each ref and alts

        """
        #Handling division by zero, when there is no ref
        for sample in self.counts:
            if sum(self.counts[sample]) <= mindepth:
                self.ratios[sample]=[math.nan]*len(self.counts[sample])
            else:
                self.ratios[sample]=[ad/sum(self.counts[sample]) for ad in self.counts[sample]]
        return self.ratios

    def scores_from_ratios(self, maxprop, minprop):
        """
        Get the scores for each alt of the variant

        maxprop (float) : Value above which a ratio is not considered a True ref
        minprop (float) : Value below which a ratio is not considered a True alt

        return (list) : Score for each alt of the variant

        """
        assert self.ratios != {}
        minprop=1-maxprop if minprop is None else minprop

        assert 0 < maxprop < 0.5 and 0.5 < minprop < 1
        for rank in range(1,len(self.alts)+1):
            ratioRank=[self.ratios[sample][rank] for sample in self.ratios]
            #No sample is above depth
            if all(np.isnan(ratio) for ratio in ratioRank):
                self.scores.append(math.nan)
                continue
            minratio=np.nanmin(ratioRank)
            maxratio=np.nanmax(ratioRank)
            #fixedVariants
            if minratio > minprop:
                self.scores.append(0)
            #Ambiguous variants
            elif minratio >= maxprop or maxratio <= minprop:
                self.scores.append(-maxratio/minratio if minratio != 0 else self.minimum)
            #True Variants
            elif minratio < maxprop and maxratio > minprop:
                self.scores.append(maxratio/minratio if minratio != 0 else self.maximum)
        return self.scores