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
        props (list) : Alt/Ref proportions of groups 1 and 2 for each alt
        miss (list) : Mssing call proportions for groups 1 and 2 for each alt
        heteroz (list) : Heterozygous call proportions for groups 1 and 2 for each alt
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
        self.counts = {}
        #AD has to be int to be considered, float not permitted (float would be set to 0)
        for i,n in enumerate(samplesRanks):
            if samples[i] in self.group1+self.group2:
                adsplit = vcfLine[n].split(":")[adRank].strip("\n").split(",")
                self.counts[samples[i]] = [int(ad) if ad.isdigit() else 0 for ad in adsplit]
                if len(self.counts[samples[i]]) < len(self.alts) + 1:
                    Config.error_print("Expected %s AD counts, only got %s. AD field was: '%s'. Set all AD counts to 0."%
                                       (len(self.alts) + 1, 
                                        len(self.counts[samples[i]]),
                                        vcfLine[n].split(":")[adRank].strip("\n")))
                    self.counts[samples[i]] = [0]*(len(self.alts) + 1)
        self.asps = {}
        self.props = []
        self.miss = []
        self.heteroz = []
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
    
    def asp_stats(self, asps):
        """
        Calculate allele population frequencies and other similar stats about samples

        Args:
        asps (list) : Alternate allele sample proportions for each sample

        return (list) : Numbers and proportions of samples
        """
        
        
        #Alternate ASP above which a loci is Alt homozygous
        minaltasp = self.config.options["minaltasp"]
        #Alternate ASP below which a loci is non Alt homozygous (assumed to be Ref homozygous)
        maxrefasp = self.config.options["maxrefasp"]
            
        lensupp = len([supptominaltasp for supptominaltasp in asps if supptominaltasp >= minaltasp])
        leninfe = len([infetomaxrefasp for infetomaxrefasp in asps if infetomaxrefasp <= maxrefasp])
        avail = lensupp+leninfe if lensupp+leninfe > 0 else math.nan
            
        missing = len([asp for asp in asps if np.isnan(asp)])
        heteroz = len(asps) - missing - avail
                       
        propsupp = lensupp/(avail + heteroz)
        propinfe = leninfe/(avail + heteroz)
        propmissing = missing/len(asps)
        propheteroz = heteroz/(avail + heteroz)
            
        maf = min(propsupp, propinfe)
        
        return [lensupp, leninfe, propsupp, propinfe, propmissing, propheteroz, maf]

    def apf_from_asps(self):
        """
        Get the mutated and reference Allele Population Frequencies using Allele Sample Proportions

        """
        assert self.asps != {}
        #Minimal and maximal population MAFs
        minmaf1group = self.config.options["minMaf1"]
        maxmaf1group = self.config.options["maxMaf1"]
        minmaf2groups = self.config.options["minMaf2"]
        maxmaf2groups = self.config.options["maxMaf2"]
        #Minimal difference of allele population frequency between groups
        minapfdiff = self.config.options["minApfDiff"]
        #Maximal proportion of missing or mixed asp in both groups
        maxmissing = self.config.options["maxMissing"]
        maxheteroz = self.config.options["maxHeteroz"]
        assert minmaf2groups <= minmaf1group <= maxmaf1group <= maxmaf2groups
        
        for rank in range(1,len(self.alts)+1):
            aspGroup1Rank = [self.asps[sample][rank] for sample in self.group1]
            aspGroup2Rank = [self.asps[sample][rank] for sample in self.group2]
            leng1supp, leng1infe, propg1supp, propg1infe, propg1missing, propg1heteroz, g1maf = self.asp_stats(aspGroup1Rank)
            leng2supp, leng2infe, propg2supp, propg2infe, propg2missing, propg2heteroz, g2maf = self.asp_stats(aspGroup2Rank)
            
            self.heteroz.append([propg1heteroz, propg2heteroz])
            self.props.append([propg1supp, propg1infe, propg2supp, propg2infe])
            self.miss.append([propg1missing, propg2missing])
            
            #Limit to number of samples in each group if 'diffsamples' too high
            #Minimal number of samples with a distinct allele in each group
            diffsamplecondition = (leng1supp >= self.diffsamples1 and leng2infe >= self.diffsamples2) or (leng1infe >= self.diffsamples1 and leng2supp >= self.diffsamples2)
            apfcondition = abs(propg1infe - propg2infe) >= minapfdiff or abs(propg1supp - propg2supp) >= minapfdiff
            
            if np.isnan(g1maf) and np.isnan(g2maf):
                mafcondition = math.nan
            else:
                maxmafcondition = maxmaf2groups >= np.nanmax([g1maf,g2maf]) >= minmaf1group
                minmafcondition = minmaf2groups <= np.nanmin([g1maf,g2maf]) <= maxmaf1group
                mafcondition = maxmafcondition and minmafcondition
            
            if (diffsamplecondition and 
                apfcondition and 
                np.nanmax([propg1heteroz, propg2heteroz]) <= maxheteroz and 
                np.nanmax([propg1missing, propg2missing]) <= maxmissing and 
                mafcondition):
                self.types.append("differential")
            elif propg1missing == 1 and propg2missing == 1:#No locus available at all (maybe some users would like to keep those too?)
                self.types.append("ambiguous")
            else:#Ambiguous alternate alleles
                self.types.append("ambiguous")
            
