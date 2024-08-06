import math
import numpy as np
from .config import Config

class Variant(object):
    """A single line of a VCF file, requiring allele depth for each sample at each allele."""
    
    def __init__(self, vcfLine, ranks, samples, samplesRanks, group1, group2, diffsamples1, diffsamples2):
        """
        Arguments:
        vcfLine (str) : Raw VCF line
        ranks (list) : Order of fields as in VCF specs
        samples (list) : Names of samples 
        samplesRanks (list) : Order of samples given in argument
        group1 (list) : Names of samples in first group of the comparison
        group2 (list) : Names of samples in second group of the comparison
        diffsamples1 (int) : Number of samples with the same allele in the first group (not the same as in second group)
        diffsamples2 (int) : Number of samples with the same allele in the second group (not the same as in first group)

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
        types (list) : Type of each alternate allele among differential and ambiguous
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
        self.diffsamples1 = diffsamples1
        self.diffsamples2 = diffsamples2
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

    def apf_from_asps(self):
        """
        Get the mutated and reference Allele Population Frequencies using Allele Sample Proportions

        """
        assert self.asps != {}
        #Value below which an asp is not considered a True mutation
        minaltasp = self.config.options["minaltasp"]
        #Value above which an asp is not considered a True reference
        maxrefasp = self.config.options["maxrefasp"]
        #Minimal difference of allele population frequency between groups
        minapfdiff = self.config.options["minApfDiff"]
        #Maximal proportion of missing or mixed asp in both groups
        maxmissing = self.config.options["maxMissing"]
        #Minimal and maximal population MAFs
        minmaf1group = self.config.options["minMaf1"]
        maxmaf1group = self.config.options["maxMaf1"]
        minmaf2groups = self.config.options["minMaf2"]
        maxmaf2groups = self.config.options["maxMaf2"]
        assert minmaf2groups <= minmaf1group <= maxmaf1group <= maxmaf2groups
        
        for rank in range(1,len(self.alts)+1):
            aspGroup1Rank = [self.asps[sample][rank] for sample in self.group1]
            aspGroup2Rank = [self.asps[sample][rank] for sample in self.group2]
            
            leng1supp = len([supptominaltasp for supptominaltasp in aspGroup1Rank if supptominaltasp >= minaltasp])
            leng1infe = len([infetomaxrefasp for infetomaxrefasp in aspGroup1Rank if infetomaxrefasp <= maxrefasp])
            leng2supp = len([supptominaltasp for supptominaltasp in aspGroup2Rank if supptominaltasp >= minaltasp])
            leng2infe = len([infetomaxrefasp for infetomaxrefasp in aspGroup2Rank if infetomaxrefasp <= maxrefasp])
            propg1supp = leng1supp/len(aspGroup1Rank)
            propg1infe = leng1infe/len(aspGroup1Rank)
            propg2supp = leng2supp/len(aspGroup2Rank)
            propg2infe = leng2infe/len(aspGroup2Rank)
            
            g1maf = min(leng1supp, leng1infe)/(leng1supp+leng1infe) if leng1supp+leng1infe > 0 else math.nan
            g2maf = min(leng2supp, leng2infe)/(leng2supp+leng2infe) if leng2supp+leng2infe > 0 else math.nan
            
            self.props.append([round(100*propg1supp), round(100*propg1infe), round(100*propg2supp), round(100*propg2infe)])

            missingg1 = 1 - propg1infe - propg1supp
            missingg2 = 1 - propg2infe - propg2supp
            apfdiffinfe = abs(propg1infe - propg2infe)
            apfdiffsupp = abs(propg1supp - propg2supp)
            
            #Limit to number of samples in each group if 'diffsamples' too high
            #Minimal number of samples with a distinct allele in each group
            diffsamplecondition = leng1supp >= self.diffsamples1 and leng2infe >= self.diffsamples2 or leng1infe >= self.diffsamples1 and leng2supp >= self.diffsamples2
            apfcondition = apfdiffinfe >= minapfdiff or apfdiffsupp >= minapfdiff
            if np.isnan(g1maf) and np.isnan(g2maf):
                mafcondition = math.nan
            else:
                maxmafcondition = maxmaf2groups >= np.nanmax([g1maf,g2maf]) >= minmaf1group
                minmafcondition = minmaf2groups <= np.nanmin([g1maf,g2maf]) <= maxmaf1group
                mafcondition = maxmafcondition and minmafcondition
            
            if diffsamplecondition and apfcondition and max(missingg1, missingg2) <= maxmissing and mafcondition:
                self.types.append("differential")
            else:#Ambiguous alternate alleles
                self.types.append("ambiguous")
            
