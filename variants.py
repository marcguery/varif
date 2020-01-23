from variant import Variant

class Variants(object):
    def __init__(self):
        self.variants={}
        self.vcfHeaderExpected=[
            "#CHROM", "POS", "ID",
            "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT"]
        self.vcfHeaderSorted={colname:i for i,colname in enumerate(self.vcfHeaderExpected)}
        self.samples=[]

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
        try:
            assert set(self.vcfHeaderExpected).difference(set(header))==set()
            self.vcfHeaderSorted={colname.strip("\n"):i for i,colname in enumerate(header)}
            n=len(header)-1
            colname=header[n].rstrip("\n")
            while colname != self.vcfHeaderExpected[-1] and n >= 0:
                self.samples.append(colname)
                n-=1
                colname=header[n]
        except Exception:
            print("Bad VCF header")

    def load_variants_from_VCF(self, vcf):
        vcffile=open(vcf, 'r').readlines()
        n=0
        line=vcffile[n]
        while line.startswith('##') and n < len(vcffile):
            line=vcffile[n]
            n+=1
        headerline=vcffile[n-1]
        self.check_vcf_header(headerline.split("\t"))
        samplesRanks=[self.vcfHeaderSorted[sample] for sample in self.samples]
        while n < len(vcffile):            
            variant=Variant(vcffile[n], self.ranks, self.samples, samplesRanks)
            variant.calculate_ratios()
            variant.scores_from_ratios()

            identifier=variant.chromosome+":"+variant.position
            num=0
            while identifier in self.variants:
                num+=1
                identifier=identifier.split(".")[0]+"."+str(num)
            self.variants[identifier]={
                "scores":variant.scores, "category":variant.category,
                "ref":variant.ref, "alts":variant.alts, 
                "ratios":variant.ratios,
                "aaRef":{}, "aaAlts":{}}
            n+=1
