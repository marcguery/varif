from sys import stderr

class Vcf(object):
    """Store data loaded from a VCF file"""

    def __init__(self):
        """
        vcfHeaderExpected (list) : Order expected of the fields in the VCF
        vcfHeaderSorted (dict) : Key (Name of field), value (order of field)
        vcffile (list) : The list of lines in the input VCF file
        headerlinenumber : The line number of the VCF file corresponding to the end of the header
        ranks (list) : Rank of each field in VCF
        samples (list) : Names of samples

        """
        self.vcfHeaderExpected=[
            "#CHROM", "POS", "ID",
            "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT"]
        self.vcfHeaderSorted={colname:i for i,colname in enumerate(self.vcfHeaderExpected)}
        self.vcffile=[]
        self.headerlinenumber=1        

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
        
        self.samples=[]
        

    def check_vcf_header(self, header):
        """
        Check if VCF format is expected and store sample names from the end

        header (list) : Name of fields in the actual VCF
        """
        try:
            assert set(self.vcfHeaderExpected).difference(set(header))==set(), "Bad VCF header"
            #Reset vcfHeaderSorted with actual header
            self.vcfHeaderSorted={colname.strip("\n"):i for i,colname in enumerate(header)}
            #Assuming that samples are at the end
            n=len(header)-1
            colname=header[n].rstrip("\n")
            while colname != self.vcfHeaderExpected[-1] and n >= 0:
                self.samples.append(colname)
                n-=1
                colname=header[n]
        except Exception as err:
            print(err, file = stderr)
            raise(err)
    
    def read_vcf(self, vcf):
        """
        Reads a VCF file

        vcf (str) : File path of VCF file
        """
        with open(vcf, 'r') as f:
            self.vcffile=f.readlines()
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