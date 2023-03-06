from sys import stderr
from .variant import Variant

class Variants(object):
    """All variants found from a VCF file."""
    def __init__(self, vcf):
        """
        vcf (Vcf) : VCF data loaded initially by varif 

        variants (dict) : Key (ID of the variant), value (Info on the variant)
        vcf (Vcf) : VCF data loaded initially by varif
        group1 (list) : Sample names in the first group to be compared 
        group2 (list) : Sample names in the second group to be compared 
        samples (list) : Sorted names of samples from both groups

        """
        self.variants={}
        self.vcf=vcf
        self.group1 = []
        self.group2 = []
    
    @property
    def samples(self):
        return sorted(list(set(self.group1+self.group2)))
    
    def check_samples(self):
        """
        Check if the sample names from the comparison groups are present in the VCF file

        return (list) : List of unique samples in group1 and group2

        """
        errors=[]
        warnings=[]
        groups=[self.group1, self.group2]
        if len(set(self.samples) - set(self.vcf.samples)) > 0:
            for index, group in enumerate(groups):
                groupmissing=set(group) - set(self.vcf.samples)
                if len(groupmissing) == len(group):
                    errors.append("No sample (%s) from one comparison group are present in the VCF file"%(", ".join(g for g in group)))
                elif len(groupmissing) > 0:
                    warnings.append("Some samples (%s) from one comparison group are absent in the VCF file"%(", ".join(g for g in groupmissing)))
                    groups[index] = list(set(groups[index])-set(groupmissing))

        if len(warnings) >  0:
            print("\n".join(warning for warning in warnings), file = stderr)
        if len(errors) > 0:
            raise NameError("\n".join(error for error in errors))
        return(groups)

    def process_variants(self, group1 = [], group2 = []):
        """
        Store variant found at each line of VCF

        group1 (list) : Sample names in the first group to be compared 
        group2 (list) : Sample names in the second group to be compared

        """
        self.group1=group1 if len(group1) > 0 else self.vcf.samples
        self.group2=group2 if len(group2) > 0 else self.vcf.samples
        self.group1, self.group2 = self.check_samples()
        samplesRanks=[self.vcf.vcfHeaderSorted[sample] for sample in self.vcf.samples]
        n = self.vcf.headerlinenumber
        while n < len(self.vcf.vcffile):
            variant=Variant(self.vcf.vcffile[n], self.vcf.ranks, self.vcf.samples, samplesRanks, self.group1, self.group2)
            variant.calculate_asps()
            variant.app_from_asps()

            #Generate unique ID for each variant
            identifier=variant.chromosome+":"+variant.position
            vcfline=n+1
            #For same Chromosome/position variants 
            # which are on same line or on different lines
            num=0
            while identifier in self.variants:
                num+=1
                identifier=identifier.split(".")[0]+"."+str(num)
            self.variants[identifier]={
                "chromosome":variant.chromosome, "position":int(variant.position),
                "props":variant.props, "types":variant.types,
                "categories":variant.categories,
                "ref":variant.ref, "alts":variant.alts,
                "refwindow":variant.refwindow,
                "asps":variant.asps,
                "features":[],
                "aaPosRef":{},
                "aaPosAlts":{},
                "cdsRef":{}, "cdsAlts":{},
                "aaRef":{}, "aaAlts":{},
                "vcfline":vcfline,
                "log":variant.log}
            n+=1
        
