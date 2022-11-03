import math
import numpy as np
from .config import Config

class Variant(object):
    """A parsed line of a VCF, requiring AD for each sample."""
    
    def __init__(self, vcfLine, ranks, samples, samplesRanks, group1, group2):
        """
        Arguments:
        vcfLine (str) : Raw VCF line
        ranks (list) : Order of fields as in VCF specs
        samples (list) : Names of samples 
        samplesRanks (list) : Order of samples given in argument
        group1 (list) : Names of samples in first group of the comparison
        group2 (list) : Names of samples in second group of the comparison

        config (Config) : Existing configuration loaded initially by varif
        chromosome (str) : Name of chromosome
        position (str) : Position of variant
        ref (str) : Reference sequence of variant starting at position
        refwindow (list) : Bases before and after the reference sequence indicated in the VCF file
        alts (list) : Alternate sequences at position
        group1 (list) : Names of samples in first group of the comparison
        group2 (list) : Names of samples in second group of the comparison
        counts (dict) : Key (sample name), value (AD count for ref and alts)
        ratios (dict) : Key (sample name), value (ratio of AD for ref and alts)
        props (list) : True var and true ref prct of both groups for each alt
        types (list) : Type of each variant among differential, fixed and ambiguous
        category (str) : Type of variant : 'SNP' or 'INDEL' and annotation if available
        log (list) : Log of the full variant (VCF line) and of each alt

        """
        assert (len(group1) == 0 and len(group2)) == 0 or (len(group1) > 0 and len(group2) > 0)
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
        self.group1=group1 if len(group1) > 0 else samples
        self.group2=group2 if len(group2) > 0 else samples
        self.counts={samples[i]:[int(ad) for ad in vcfLine[n].split(":")[adRank].strip("\n").split(",")] for i,n in enumerate(samplesRanks) if samples[i] in self.group1+self.group2}
        self.ratios={}
        self.props=[]
        self.types=[]
        self.category=""   
        self.log=["",[]]    
    
    def calculate_ratios(self):
        """
        Get the ratio for each alternate count of the variant

        return (dict) : Ratios for each ref and alts

        """
        try:
            mindepth=int(self.config.options["mindepth"])
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
            ratioGroup1Rank=[self.ratios[sample][rank] for sample in self.group1]
            ratioGroup2Rank=[self.ratios[sample][rank] for sample in self.group2]
            #No sample is above depth
            if all(np.isnan(ratio) for ratio in ratioGroup1Rank):
                ming1ratio=math.nan
                maxg1ratio=math.nan
            else:
                ming1ratio=np.nanmin(ratioGroup1Rank)
                maxg1ratio=np.nanmax(ratioGroup1Rank)
            if all(np.isnan(ratio) for ratio in ratioGroup2Rank):
                ming2ratio=math.nan
                maxg2ratio=math.nan
            else:
                ming2ratio=np.nanmin(ratioGroup2Rank)
                maxg2ratio=np.nanmax(ratioGroup2Rank)
            
            propg1supp=len([supptominaaf for supptominaaf in ratioGroup1Rank if supptominaaf>=minaaf])/len(ratioGroup1Rank)
            propg1infe=len([infetomaxraf for infetomaxraf in ratioGroup1Rank if infetomaxraf<=maxraf])/len(ratioGroup1Rank)
            propg2supp=len([supptominaaf for supptominaaf in ratioGroup2Rank if supptominaaf>=minaaf])/len(ratioGroup2Rank)
            propg2infe=len([infetomaxraf for infetomaxraf in ratioGroup2Rank if infetomaxraf<=maxraf])/len(ratioGroup2Rank)
            
            self.props.append([round(100*propg1supp), round(100*propg1infe), round(100*propg2supp), round(100*propg2infe)])

            missingg1 = 1 - propg1infe - propg1supp
            missingg2 = 1 - propg2infe - propg2supp
            #fixed alternate or reference
            if ming1ratio >= minaaf and ming2ratio >= minaaf or (maxg1ratio <= maxraf and maxg2ratio <= maxraf):
                if max(missingg1, missingg2) <= maxmissing:
                    self.types.append("fixed")
                else:
                    self.types.append("ambiguous")
            #True Alt
            elif ming1ratio <= maxraf and maxg2ratio >= minaaf or (ming2ratio <= maxraf and maxg1ratio >= minaaf):
                propvariantdiff = min(propg1supp,propg2supp) / max(propg1supp,propg2supp)
                if max(missingg1, missingg2) <= maxmissing and max(propg1supp,propg2supp) >= minvariants and propvariantdiff <= ratiovariantdiff:
                    self.types.append("differential")
                else:
                    self.types.append("ambiguous")
            #Ambiguous variants
            else:
                self.types.append("ambiguous")
