import math
import numpy as np
from .config import Config

class Variant(object):
    """A parsed line of a VCF, requiring AD for each sample."""
    def __init__(self, vcfLine, ranks, samples, samplesRanks, refSamples):
        """
        Arguments:
        vcfLine (str) : Raw VCF line
        ranks (list) : Order of fields as in VCF specs
        samples (list) : Name of samples 
        samplesRanks (list) : Order of samples given in argument
        refSamples (list) : Names of grouped samples

        chromosome (str) : Name of chromosome
        position (str) : Position of variant
        ref (str) : Reference sequence of variant starting at position
        alts (list) : Alternate sequences at position
        counts (dict) : Key (sample name), value (AD count for ref and alts)
        ratios (dict) : Key (sample name), value (ratio of AD for ref and alts)
        props (list) : True var and true ref prct of both groups for each alt
        types (list) : Type of each variant among differential, fixed and ambiguous
        category (str) : Type of variant : 'SNP' or 'INDEL' and annotation if available
        log (list) : Log of the full variant (VCF line) and of each alt

        """
        self.config=Config
        vcfLine=vcfLine.split("\t")
        self.chromosome=vcfLine[ranks[0]]
        self.position=vcfLine[ranks[1]]
            
        self.ref=vcfLine[ranks[3]]
        self.refwindow=["",""]
        self.alts=vcfLine[ranks[4]].split(",")
        formatSplitted=vcfLine[ranks[8]].split(":")
        for i in range(len(formatSplitted)):
            if formatSplitted[i] == "AD":
                adRank=i
                break
        self.counts={samples[i]:[int(ad) for ad in vcfLine[n].split(":")[adRank].strip("\n").split(",")] for i,n in enumerate(samplesRanks)}
        self.refSamples=refSamples if len(refSamples) > 0 else samples
        self.altSamples=list(set(samples)-set(self.refSamples)) if len(refSamples) > 0 else samples
        self.ratios={}
        self.props=[]
        self.types=[]
        self.category=""   
        self.log=["",[]]    
    
    def calculate_ratios(self, mindepth):
        """
        Get the ratio for each alternate count of the variant

        mindepth (int) : Minmal number of reads (REF+ALTs) to
        calculate allele frequencies, must be integer >= 1

        return (dict) : Ratios for each ref and alts

        """
        try:
            mindepth=int(mindepth)
        except ValueError:
            raise ValueError("Argument 'mindepth' should be an integer, not '%s'"%mindepth)
        if mindepth < 1:
            raise ValueError("Depth should be above 0")
        #Handling division by zero, when there is no ref
        for sample in self.counts:
            if sum(self.counts[sample]) < mindepth:
                self.ratios[sample]=[math.nan]*len(self.counts[sample])
            else:
                self.ratios[sample]=[round(ad/sum(self.counts[sample]), 2) for ad in self.counts[sample]]
        return self.ratios

    def props_from_ratios(self):
        """
        Get the true var and true ref prct for each alt of the variant

        """
        assert self.ratios != {}
        #Value below which a ratio is not considered a True alt
        minaaf=self.config.options["minaaf"]
        #Value above which a ratio is not considered a True ref
        maxraf=self.config.options["maxraf"]
        minaaf=1-maxraf if minaaf is None else minaaf

        #Maximal proportion of missing (NA or mixed AF) in both groups
        maxmissing = self.config.options["maxMissing"]
        #Minimal proportion of samples that are called true alt (differential groups only)
        minvariants = self.config.options["minVariants"]
        #Maximal ratio of the least number of alt samples over the other group alt samples (differential groups only)
        ratiovariantdiff = self.config.options["maxSimilarity"]
        assert 0 <= maxraf <= 1 and 0 <= minaaf <= 1, "Alternate AF (input:%s) and Reference AF (input:%s) must be proportions"%(minaaf, maxraf)
        assert maxraf <= minaaf, "Reference AF should be <= to Alternate AF"
        assert 0 <= maxmissing <= 1 and 0 <= minvariants <= 1 and 0 <= ratiovariantdiff <= 1, "Group filtering options must be proportions"
        
        for rank in range(1,len(self.alts)+1):
            ratioRefRank=[self.ratios[sample][rank] for sample in self.refSamples]
            ratioAltRank=[self.ratios[sample][rank] for sample in self.altSamples]
            #No sample is above depth
            if all(np.isnan(ratio) for ratio in ratioRefRank):
                minrefratio=math.nan
                maxrefratio=math.nan
            else:
                minrefratio=np.nanmin(ratioRefRank)
                maxrefratio=np.nanmax(ratioRefRank)
            if all(np.isnan(ratio) for ratio in ratioAltRank):
                minaltratio=math.nan
                maxaltratio=math.nan
            else:
                minaltratio=np.nanmin(ratioAltRank)
                maxaltratio=np.nanmax(ratioAltRank)

            numaltsupp=len([supptominaaf for supptominaaf in ratioAltRank if supptominaaf>=minaaf])/len(ratioAltRank)
            numaltinfe=len([infetomaxraf for infetomaxraf in ratioAltRank if infetomaxraf<=maxraf])/len(ratioAltRank)
            numrefsupp=len([supptominaaf for supptominaaf in ratioRefRank if supptominaaf>=minaaf])/len(ratioRefRank)
            numrefinfe=len([infetomaxraf for infetomaxraf in ratioRefRank if infetomaxraf<=maxraf])/len(ratioRefRank)
            self.props.append([round(100*numaltsupp), round(100*numaltinfe), round(100*numrefsupp), round(100*numrefinfe)])

            missingref = 1 - numrefinfe - numrefsupp
            missingalt = 1 - numaltinfe - numaltsupp
            #fixed alternate or reference
            if minrefratio >= minaaf and minaltratio >= minaaf or (maxrefratio <= maxraf and maxaltratio <= maxraf):
                if max(missingref, missingalt) <= maxmissing:
                    self.types.append("fixed")
                else:
                    self.types.append("ambiguous")
            #True Alt
            elif minrefratio <= maxraf and maxaltratio >= minaaf or (minaltratio <= maxraf and maxrefratio >= minaaf):
                propvariantdiff = min(numrefsupp,numaltsupp) / max(numrefsupp,numaltsupp)
                if max(missingref, missingalt) <= maxmissing and max(numrefsupp,numaltsupp) >= minvariants and propvariantdiff <= ratiovariantdiff:
                    self.types.append("differential")
                else:
                    self.types.append("ambiguous")
            #Ambiguous variants
            else:
                self.types.append("ambiguous")
