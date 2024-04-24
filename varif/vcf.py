from sys import stderr
import os

class Vcf(object):
    """Store data loaded from a VCF file"""

    vcfpath=""
    vcfHeaderExpected=[
            "#CHROM", "POS", "ID",
            "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT"]
    vcfHeaderSorted = {}
    headerlinenumber = -1
    totlines = -1
    intervals = []
    ranks = []
    samples = []

    def __init__(self):
        """
        vcfpath (str) : Path to the VCf file
        vcfHeaderExpected (list) : Order expected of the fields in the VCF
        vcfHeaderSorted (dict) : Key (Name of field), value (order of field)
        headerlinenumber : The line number of the VCF file corresponding to the end of the header
        totlines : The total number of lines in the VCf file
        intervals (list) : Ranges of VCF lines to use for each chunk
        ranks (list) : Rank of each field in VCF
        samples (list) : Names of samples
        vcffile (list) : The list of lines in the input VCF file

        """
        self.vcffile=[]
    
    @classmethod
    def check_vcf_header(self, header):
        """
        Check if VCF format is expected and store VCF fields and sample names

        header (list) : Name of fields in the actual VCF
        """
        try:
            assert set(self.vcfHeaderExpected).difference(set(header))==set(), "Bad VCF header"
            self.vcfHeaderSorted={colname:i for i,colname in enumerate(self.vcfHeaderExpected)}
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
    
    @classmethod
    def _read_vcf_header(self, vcf):
        """
        Reads a VCF header and store its location

        vcf (str) : File path of VCF file

        return (str) : The lines of the header

        """
        self.vcfpath = vcf
        vcffile = []
        with open(self.vcfpath, 'r') as f:
            line = f.readline()
            vcffile.append(line)
            n=1
            while line.startswith('##'):
                line = f.readline()
                vcffile.append(line)
                n+=1
        self.headerlinenumber=n
        headerline=line
        #Check VCF formatting
        self.check_vcf_header(headerline.split("\t"))
        return vcffile
    
    def read_vcf_header(self, vcf):
        """
        Save a VCF header

        vcf (str) : File path of VCF file

        """
        self.vcffile = self._read_vcf_header(vcf)
    
    @classmethod
    def count_lines(self):
        """Count the number of lines in a VCF file."""
        self.totlines = sum(1 for _ in open(self.vcfpath))
    
    @classmethod
    def define_intervals(self, chunks):
        """
        Determine the ranges of lines in the VCF file matching each chunk

        chunks (int) : The number of chunks to divide the VCF file into

        """
        optimalchunks = chunks if (self.totlines-self.headerlinenumber)%chunks == 0 else chunks-1
        chunk = (self.totlines-self.headerlinenumber)//optimalchunks
        self.intervals = [[self.headerlinenumber+part*chunk,self.headerlinenumber+(part+1)*chunk] for part in range(optimalchunks)]
        if (self.totlines-self.headerlinenumber)%chunks > 0:
            self.intervals.append([self.headerlinenumber+optimalchunks*chunk, self.headerlinenumber+optimalchunks*chunk + (self.totlines-self.headerlinenumber)%optimalchunks])
        averagechunksize = int(sum([int[1]-int[0] for int in self.intervals])/len(self.intervals))
        if  averagechunksize < 450:
            print("Low average chunk size detected of %s, consider increasing the chunk size upper limit or reducing the number of cores."%averagechunksize)
        
class Vcfdata(Vcf):
    """Store variant data loaded from a VCF file"""

    def __init__(self, chunknumber = 1):
        """
        vcffile (list) : The content of the VCF file
        chunknumber (int) : The rank of the chunk of variants to use
        interval (list) : The line numbers where the variants are located

        """
        #The VCF file has to be read several times, it should not be spawned by a subprocess
        if not os.path.isfile(self.vcfpath):
            raise FileNotFoundError("VCF file '%s' should be a regular file."%self.vcfpath)
        
        self.vcffile=[]
        self.chunknumber = chunknumber
    
    @property
    def interval(self):
        return [self.intervals[self.chunknumber-1][0],self.intervals[self.chunknumber-1][1]]
    
    def read_vcf(self):
        """Reads variants from a VCF file"""
        firstline=self.interval[0]
        lastline=self.interval[1]
        self.vcffile=[]
        currentline=0
        with open(self.vcfpath, 'r') as f:
            while currentline < lastline:
                currentline+=1
                if currentline > max(self.headerlinenumber, firstline):
                    self.vcffile.append(f.readline())
                else:
                    f.readline()
