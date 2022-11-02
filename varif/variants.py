from sys import stderr
from .config import Config
from .variant import Variant

class Variants(object):
    """A bunch of variants taken from a VCF file."""
    def __init__(self):
        """
        variants (dict) : Key (ID of the variant), value (Info on the variant)
        vcfHeaderExpected (list) : Order expected of the fields in the VCF
        vcfHeaderSorted (dict) : Key (Name of field), value (order of field)
        vcffile (list) : The list of lines in the input VCF file
        headerlinenumber : The line number of the VCF file corresponding to the end of the header
        samples (list) : Name of samples
        refRanks (list) : Rank of each grouped sample (starting at 1)
        refSamples (list) : Name of each grouped sample
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
        refRanksraw=Config.options["refRanks"]
        self.refRanks=[]
        for arg in refRanksraw.split(","):
            #In the case there is no grouped samples registered:
            if arg == "" and len(refRanksraw.split(",")) == 1:
                break
            try:
                self.refRanks.append(int(arg.strip("")))
            except Exception as err:
                print("Bad format of the grouped sample numbers", file = stderr)
                print(err, file = stderr)
        self.refSamples=[]

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
            assert set(self.vcfHeaderExpected).difference(set(header))==set(), "Bad VCF header"
            if self.refRanks != []:
                assert min(self.refRanks) >= 1 and max(self.refRanks) <= len(header)-len(self.vcfHeaderExpected), "Wrong sample ID"
            #Reset vcfHeaderSorted with actual header
            self.vcfHeaderSorted={colname.strip("\n"):i for i,colname in enumerate(header)}
            #Assuming that samples are at the end
            n=len(header)-1
            colname=header[n].rstrip("\n")
            while colname != self.vcfHeaderExpected[-1] and n >= 0:
                self.samples.append(colname)
                if n - (len(self.vcfHeaderExpected)-1) in self.refRanks:
                    self.refSamples.append(colname)
                n-=1
                colname=header[n]
        except Exception as err:
            print(err, file = stderr)
            raise(err)

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
        if len(self.refSamples) > 0:
            print("Grouped samples are : %s"%(", ".join(self.refSamples)))
        #Storing variants
        while n < len(self.vcffile):
            variant=Variant(self.vcffile[n], self.ranks, self.samples, samplesRanks, self.refSamples)
            variant.calculate_ratios(self.config.options["mindepth"])
            variant.props_from_ratios()

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
                "category":variant.category,
                "ref":variant.ref, "alts":variant.alts,
                "refwindow":variant.refwindow,
                "ratios":variant.ratios,
                "features":[],
                "aaPosRef":{},
                "aaPosAlts":{},
                "cdsRef":{}, "cdsAlts":{},
                "aaRef":{}, "aaAlts":{},
                "vcfline":vcfline,
                "log":variant.log}
            n+=1
