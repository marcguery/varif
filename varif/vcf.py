from sys import stderr
import os

class Vcf(object):
    """Store data loaded from a VCF file"""

    vcfpath = ""
    vcfheader = []
    vcfHeaderExpected = [
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
        vcfpath (str) : Path to the VCF file
        vcfheader (list) : List of lines in the VCF herader (up tp "#CHROM")
        vcfHeaderExpected (list) : Order expected of the fields in the VCF
        vcfHeaderSorted (dict) : Key (Name of field), value (order of field)
        headerlinenumber : The line number of the VCF file corresponding to the end of the header
        totlines : The total number of lines in the VCF file
        intervals (list) : Ranges of VCF lines to use for each chunk
        ranks (list) : Rank of each field in VCF
        samples (list) : Names of samples
        vcffile (list) : List of lines after the header (after "#CHROM") in the VCF

        """
        self.vcffile = []
    
    @classmethod
    def check_vcf_header(cls, header):
        """
        Check if VCF format is expected and store VCF fields and sample names

        header (list) : Fields of the VCF
        
        """
        try:
            assert set(cls.vcfHeaderExpected).difference(set(header)) == set(), "Bad VCF header"
            cls.vcfHeaderSorted = {colname:i for i,colname in enumerate(cls.vcfHeaderExpected)}
            #A very complicated stuff where a simple order to follow would have been enough
            chromRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[0]]
            posRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[1]]
            idRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[2]]
            refRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[3]]
            altsRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[4]]
            qualRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[5]]
            filterRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[6]]
            infoRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[7]]
            formatRank = cls.vcfHeaderSorted[cls.vcfHeaderExpected[8]]
            cls.ranks = [
            chromRank, posRank, idRank,
            refRank, altsRank, qualRank,
            filterRank, infoRank, formatRank]
            #Reset vcfHeaderSorted with actual header
            cls.vcfHeaderSorted = {colname.strip("\n"):i for i,colname in enumerate(header)}
            #Assuming that samples are at the end
            n = len(header)-1
            colname = header[n].rstrip("\n")
            while colname != cls.vcfHeaderExpected[-1] and n >= 0:
                cls.samples.append(colname)
                n -= 1
                colname = header[n]
        except Exception as err:
            print(err, file = stderr)
            raise(err)
    
    @classmethod
    def read_vcf_header(cls, vcf):
        """
        Reads a VCF header and store its location

        vcf (str) : File path of VCF file

        """
        cls.vcfpath = vcf
        cls.vcfheader = []
        with open(cls.vcfpath, 'r') as f:
            n = 0
            line = f.readline()
            while line.startswith('##'):
                n += 1
                cls.vcfheader.append(line)
                line = f.readline()
            if line.startswith(cls.vcfHeaderExpected[0]):
                n += 1
                cls.vcfheader.append(line)
        cls.headerlinenumber = n
        #Check VCF formatting
        cls.check_vcf_header(cls.vcfheader[-1].split("\t"))
    
    @classmethod
    def count_lines(cls):
        """Count the number of lines in a VCF file."""
        cls.totlines = sum(1 for _ in open(cls.vcfpath))
    
    @classmethod
    def define_intervals(cls, chunks):
        """
        Determine the evenly distributed ranges of lines in the VCF file matching each chunk

        chunks (int) : The number of chunks to divide the VCF file into
        
        return (list) : The unique sizes of all intervals

        """
        assert (cls.totlines-cls.headerlinenumber) >= chunks >= 1
        
        #To make the intervals as even as possible, they are formed using the formula:
        # n = i * (x+1) + j * x
        # with:
        # n = total number of variants,
        # x = n//(# of chunks),
        # i and j the numbers of intervals of size x or x+1 so that the sum equals n
        
        base_interval = (cls.totlines - cls.headerlinenumber)//chunks #corresponds to 'x'
        last_interval = cls.totlines - cls.headerlinenumber - (base_interval * (chunks - 1))
        extra = last_interval - base_interval #remainder is equal to 'i'
        majorshift = cls.headerlinenumber #intervals must consider the VCF header as well

        cls.intervals = [[majorshift+part*(base_interval+1),majorshift+(part+1)*(base_interval+1)] for part in range(extra)]
        cls.intervals += [[majorshift+extra+part*(base_interval),majorshift+extra+(part+1)*(base_interval)] for part in range(extra,chunks)]
        
        uniquechunksizes = [size_inter for size_inter in set([interval[1]-interval[0] for interval in cls.intervals])]
        if  min(uniquechunksizes) < 100:
            log = ("Low chunk size of %s detected,"
                    " consider increasing the chunk size upper limit or reducing the number of cores.")%min(uniquechunksizes)
            print(log)
        
        return uniquechunksizes
    
    def read_vcf(self):
        """Reads the whole VCF file without the header"""
        if len(self.vcfheader) == 0:
            raise Exception("You should read the VCF header first.")
        self.vcffile = []
        currentline = 0
        with open(self.vcfpath, 'r') as f:
            for line in f:
                currentline += 1
                if currentline > self.headerlinenumber:
                    self.vcffile.append(line)
            

class Vcfchunk(Vcf):
    """Store variant chunks from a VCF file"""

    def __init__(self, chunknumber = 1):
        """
        chunknumber (int) : The rank of the chunk of variants to use
        vcfchunk (list) : The list of lines in the VCF chunk
        interval (list) : The numbers of the first and last lines where the variants are located

        """
        #The VCF file has to be read several times, it should not be spawned by a subprocess
        if not os.path.isfile(self.vcfpath):
            raise FileNotFoundError("VCF file '%s' should be a regular file."%self.vcfpath)
        
        self.chunknumber = chunknumber
        self.vcfchunk = []
    
    @property
    def interval(self):
        return [self.intervals[self.chunknumber-1][0],self.intervals[self.chunknumber-1][1]]
    
    def read_vcf(self):
        """Reads variants from a chunk of a VCF file"""
        firstline = max(self.headerlinenumber, self.interval[0])
        lastline = self.interval[1]
        self.vcfchunk = []
        currentline = 0
        with open(self.vcfpath, 'r') as f:
            while currentline < lastline:
                currentline += 1
                if currentline > firstline:
                    self.vcfchunk.append(f.readline())
                else:
                    f.readline()
