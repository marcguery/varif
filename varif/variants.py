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
        self.variants = {}
        self.vcf = vcf
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
        errors = []
        warnings = []
        groups = [self.group1, self.group2]
        if len(set(self.samples) - set(self.vcf.samples)) > 0:
            for index, group in enumerate(groups):
                groupmissing = set(group) - set(self.vcf.samples)
                if len(groupmissing) == len(group):
                    errors.append("No sample from one comparison group (%s) are present in the VCF file '%s'." % (
                        ", ".join(g for g in group), self.vcf.vcfpath))
                    errors.append("Check your PED file or your VCF header.")
                elif len(groupmissing) > 0:
                    warnings.append("Some samples from one comparison group (%s) are absent in the VCF file '%s'." % (
                        ", ".join(g for g in groupmissing), self.vcf.vcfpath))
                    warnings.append("Check your PED file or your VCF header.")
                    groups[index] = list(set(groups[index])-set(groupmissing))

        if len(warnings) > 0:
            print("\n".join(warning for warning in warnings), file = stderr)
        if len(errors) > 0:
            raise NameError("\n".join(error for error in errors))
        return (groups)

    def process_vcf_variant(self, vcfline):
        """
        Calculate ASPs then APPs from variants at a given VCF line

        vcfline (int) : Variants from this line number of the VCF will be processed

        """
        variant = Variant(self.vcf.vcffile[vcfline-1], self.vcf.ranks,
                          self.vcf.samples, [self.vcf.vcfHeaderSorted[sample]
                        for sample in self.vcf.samples], self.group1, self.group2)
        variant.calculate_asps()
        variant.app_from_asps()

        # Generate unique ID for each variant
        # and for same Chromosome/position variants
        # which are on different lines
        identifier = variant.chromosome+":"+variant.position+"."+str(vcfline)
        self.variants[identifier] = {
            "chromosome": variant.chromosome, "position": int(variant.position),
            "props": variant.props, "types": variant.types,
            "categories": variant.categories,
            "ref": variant.ref, "alts": variant.alts,
            "refwindow": variant.refwindow,
            "asps": variant.asps,
            "features": [],
            "aaPosRef": {},
            "aaPosAlts": {},
            "cdsRef": {}, "cdsAlts": {},
            "aaRef": {}, "aaAlts": {},
            "vcfline": vcfline,
            "log": variant.log}

    def process_variants(self, group1 = [], group2 = []):
        """
        Store variant found at each line of VCF

        group1 (list) : Sample names in the first group to be compared 
        group2 (list) : Sample names in the second group to be compared

        """
        self.group1 = group1 if len(group1) > 0 else self.vcf.samples
        self.group2 = group2 if len(group2) > 0 else self.vcf.samples
        self.group1, self.group2 = self.check_samples()
        
        n = 1
        while n <= len(self.vcf.vcffile):
            self.process_vcf_variant(n)
            n += 1
