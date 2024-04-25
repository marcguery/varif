import math
import numpy as np
from .config import Config

class Variant(object):
    """A single line of a VCF file, requiring allele depth for each sample at each allele."""
    
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
        ref (str) : Reference allele
        refwindow (list) : Bases before and after the reference sequence to be included
        alts (list) : Alternate sequences at given position
        group1 (list) : Names of samples in first group of the comparison
        group2 (list) : Names of samples in second group of the comparison
        counts (dict) : Key (sample name), value (AD of each sample for each allele including ref)
        asps (dict) : Key (sample name), value (asp of each sample for each allele including ref)
        props (list) : True var and true ref prct of both groups for each alt
        types (list) : Type of each alternate allele among differential, fixed and ambiguous
        categories (list) : Type of each alternate allele at the given coordinates : 'SNP' or 'INDEL'
        log (list) : Log of the full variant (VCF line) and of each alt

        """
        assert (len(group1) == 0 and len(group2)) == 0 or (len(group1) > 0 and len(group2) > 0)
        self.config = Config
        vcfLine = vcfLine.split("\t")
        self.chromosome = vcfLine[ranks[0]]
        self.position = vcfLine[ranks[1]]
            
        self.ref = vcfLine[ranks[3]]
        self.refwindow = ["",""]
        self.alts = vcfLine[ranks[4]].split(",")
        formatSplitted = vcfLine[ranks[8]].split(":")
        for i in range(len(formatSplitted)):
            if formatSplitted[i] == "AD":
                adRank = i
                break
        self.group1 = group1 if len(group1) > 0 else samples
        self.group2 = group2 if len(group2) > 0 else samples
        self.counts = {samples[i]:[int(ad) for ad in vcfLine[n].split(":")[adRank].strip("\n").split(",")] for i,n in enumerate(samplesRanks) if samples[i] in self.group1+self.group2}
        self.asps = {}
        self.props = []
        self.types = []
        self.categories = []   
        self.log = ["",[]]
    
    def calculate_asps(self):
        """
        Get the Allele Sample Proportion (ASP) for each sample at each allele including ref

        return (dict) : ASPs for each sample

        """
        try:
            mindepth = int(self.config.options["mindepth"])
        except ValueError:
            raise ValueError("Argument 'mindepth' should be an integer, not '%s'"%mindepth)
        if mindepth < 1:
            raise ValueError("Depth should be above 0")
        #Handling division by zero, when there is no ref
        for sample in self.counts:
            if sum(self.counts[sample]) < mindepth:
                self.asps[sample] = [math.nan]*len(self.counts[sample])
            else:
                self.asps[sample] = [round(ad/sum(self.counts[sample]), 6) for ad in self.counts[sample]]
        return self.asps

    def app_from_asps(self):
        """
        Get the mutated and reference Allele Population Proportions using Allele Sample Proportions

        """
        assert self.asps != {}
        #Value below which an asp is not considered a True mutation
        minaltasp = self.config.options["minaltasp"]
        #Value above which an asp is not considered a True reference
        maxrefasp = self.config.options["maxrefasp"]

        #Maximal proportion of missing or mixed asp in both groups
        maxmissing = self.config.options["maxMissing"]
        #Minimal proportion of samples that are have a true mutation (differential groups only)
        minmutated = self.config.options["minMutations"]
        #Maximal ratio of min(# of mutated samples)/max(# of mutated samples) (differential groups only)
        ratiomutationdiff = self.config.options["maxSimilarity"]
        
        for rank in range(1,len(self.alts)+1):
            aspGroup1Rank = [self.asps[sample][rank] for sample in self.group1]
            aspGroup2Rank = [self.asps[sample][rank] for sample in self.group2]
            #No sample is above depth
            if all(np.isnan(asp) for asp in aspGroup1Rank):
                ming1asp = math.nan
                maxg1asp = math.nan
            else:
                ming1asp = np.nanmin(aspGroup1Rank)
                maxg1asp = np.nanmax(aspGroup1Rank)
            if all(np.isnan(asp) for asp in aspGroup2Rank):
                ming2asp = math.nan
                maxg2asp = math.nan
            else:
                ming2asp = np.nanmin(aspGroup2Rank)
                maxg2asp = np.nanmax(aspGroup2Rank)
            
            propg1supp = len([supptominaltasp for supptominaltasp in aspGroup1Rank if supptominaltasp >= minaltasp])/len(aspGroup1Rank)
            propg1infe = len([infetomaxrefasp for infetomaxrefasp in aspGroup1Rank if infetomaxrefasp <= maxrefasp])/len(aspGroup1Rank)
            propg2supp = len([supptominaltasp for supptominaltasp in aspGroup2Rank if supptominaltasp >= minaltasp])/len(aspGroup2Rank)
            propg2infe = len([infetomaxrefasp for infetomaxrefasp in aspGroup2Rank if infetomaxrefasp <= maxrefasp])/len(aspGroup2Rank)
            
            self.props.append([round(100*propg1supp), round(100*propg1infe), round(100*propg2supp), round(100*propg2infe)])

            missingg1 = 1 - propg1infe - propg1supp
            missingg2 = 1 - propg2infe - propg2supp
            #fixed mutation or reference
            if ming1asp >= minaltasp and ming2asp >= minaltasp or (maxg1asp <= maxrefasp and maxg2asp <= maxrefasp):
                if max(missingg1, missingg2) <= maxmissing:
                    self.types.append("fixed")
                else:
                    self.types.append("ambiguous")
            #True mutational difference between groups
            elif ming1asp <= maxrefasp and maxg2asp >= minaltasp or (ming2asp <= maxrefasp and maxg1asp >= minaltasp):
                propmutationdiff = min(propg1supp,propg2supp) / max(propg1supp,propg2supp)
                if max(missingg1, missingg2) <= maxmissing and max(propg1supp,propg2supp) >= minmutated and propmutationdiff <= ratiomutationdiff:
                    self.types.append("differential")
                else:
                    self.types.append("ambiguous")
            #Ambiguous alternate alleles
            else:
                self.types.append("ambiguous")
