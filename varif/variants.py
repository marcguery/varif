from .config import Config
from .variant import Variant

class Variants(object):
    """A bunch if variants taken from a VCF file."""
    def __init__(self):
        """
        variants (dict) : Key (ID of the variant), value (Info on the variant)
        vcfHeaderExpected (list) : Order expected of the fields in the VCF
        vcfHeaderSorted (dict) : Key (Name of field), value (order of field)
        vcffile (list) : The list of lines in the input VCF file
        headerlinenumber : The line number of the VCF file corresponding to the end of the header
        samples (list) : Name of samples
        ranks (list) : Rank of each field in VCF

        """
        self.variants={}
        self.vcfHeaderExpected=[
            "#CHROM", "POS", "ID",
            "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT"]
        self.vcfHeaderSorted={colname:i for i,colname in enumerate(self.vcfHeaderExpected)}
        self.vcffile=[]
        self.headerlinenumber=1
        self.samples=[]
        self.config=Config

        #A very complicated stuff where a simple order to follow would have been enough
        chromRank=self.vcfHeaderSorted[self.vcfHeaderExpected[0]]
        posRank=self.vcfHeaderSorted[self.vcfHeaderExpected[1]]
        idRank=self.vcfHeaderSorted[self.vcfHeaderExpected[2]]
        refRank=self.vcfHeaderSorted[self.vcfHeaderExpected[3]]
        altsRank=self.vcfHeaderSorted[self.vcfHeaderExpected[4]]
        qualRank=self.vcfHeaderSorted[self.vcfHeaderExpected[5]]
        filterRank=self.vcfHeaderSorted[self.vcfHeaderExpected[6]]
        infoRank=self.vcfHeaderSorted[self.vcfHeaderExpected[7]]
        formatRank=self.vcfHeaderSorted[self.vcfHeaderExpected[8]]

        self.ranks=[
            chromRank, posRank, idRank,
            refRank, altsRank, qualRank,
            filterRank, infoRank, formatRank]
    def check_vcf_header(self, header):
        """
        Check if VCF format is expected and stroe sample names from the end

        header (list) : Name of fields in the actual VCF
        """
        try:
            assert set(self.vcfHeaderExpected).difference(set(header))==set()
            #Reset vcfHeaderSorted with actual header
            self.vcfHeaderSorted={colname.strip("\n"):i for i,colname in enumerate(header)}
            #Assuming that samples are at the end
            n=len(header)-1
            colname=header[n].rstrip("\n")
            while colname != self.vcfHeaderExpected[-1] and n >= 0:
                self.samples.append(colname)
                n-=1
                colname=header[n]
        except Exception:
            print("Bad VCF header")

    def load_variants_from_VCF(self, vcf):
        """
        Store variant found at each line of VCF

        vcf (str) : File path of VCF file

        """
        self.vcffile=open(vcf, 'r').readlines()
        n=0
        line=self.vcffile[n]
        #Getting rid of the header
        while line.startswith('##') and n < len(self.vcffile):
            line=self.vcffile[n]
            n+=1
        self.headerlinenumber=n
        headerline=self.vcffile[self.headerlinenumber-1]
        #Check VCF formatting
        self.check_vcf_header(headerline.split("\t"))
        samplesRanks=[self.vcfHeaderSorted[sample] for sample in self.samples]
        #Storing variants
        while n < len(self.vcffile):
            variant=Variant(self.vcffile[n], self.ranks, self.samples, samplesRanks)
            variant.calculate_ratios(self.config.options["mindepth"])
            variant.scores_from_ratios()

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
                "scores":variant.scores, "category":variant.category,
                "ref":variant.ref, "alts":variant.alts, 
                "ratios":variant.ratios, "groups":variant.groups,
                "features":[],
                "aaRef":{}, "aaAlts":{},
                "vcfline":vcfline}
            n+=1
