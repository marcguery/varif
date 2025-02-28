from .config import Config
from os import access, R_OK
from os.path import isfile

class Vcf(object):
    """Store data loaded from a VCF file"""

    def __init__(self):
        """
        vcfpath (str) : Path to the VCF file
        vcfheader (list) : List of lines in the VCF herader (up tp "#CHROM")
        vcfHeaderExpected (list) : Order expected of the fields in the VCF
        headerlinenumber : The line number of the VCF file corresponding to the end of the header
        totlines : The total number of lines in the VCF file
        intervals (list) : Ranges of VCF lines to use for each chunk
        samples (list) : Names of samples
        vcfbody (list) : List of lines after the header (after "#CHROM") in the VCF

        """
        self.reset()
        
    def reset(self):
        """Reset attributes"""
        self.vcfpath = ""
        self.vcfheader = []
        self.vcfHeaderExpected = [
                "#CHROM", "POS", "ID",
                "REF", "ALT", "QUAL", "FILTER",
                "INFO", "FORMAT"]
        self.headerlinenumber = -1
        self.totlines = -1
        self.intervals = []
        self.samples = []
        self.vcfbody = []
        
    def readable_file(self):
        """
        Checks if the VCF file can be read

        return (bool) : If file can be read or not
        """
        return isfile(self.vcfpath) and access(self.vcfpath, R_OK)
    
    def check_vcf_header(self, header):
        """
        Check if VCF format is expected and store VCF fields and sample names

        header (list) : Fields of the VCF
        
        """
        if not set(self.vcfHeaderExpected).difference(set(header)) == set():
            Config.error_print("Unexpected VCF header")
            raise ValueError
        
        for n, colname in enumerate(self.vcfHeaderExpected):
            if header[n].rstrip("\n") != colname:
                Config.error_print("Wrong ordering of VCF header, offending column #%s named '%s', expected '%s'"%(
                    n+1, header[n], colname
                ))
                raise ValueError
        
        if len(header) == len(self.vcfHeaderExpected):
            Config.error_print("No sample found in the VCF")
            raise ValueError
        n = 9
        while n < len(header):
            self.samples.append(header[n].rstrip("\n"))
            n += 1
        Config.verbose_print("Found %s samples."%(len(self.samples)))
    
    def read_vcf(self, vcf, headerOnly = False):
        """
        Reads a VCF header and store its location

        vcf (str) : File path of VCF file
        headerOnly (bool) : Reads the header only

        """
        self.reset()
        self.vcfpath = vcf
        if not self.readable_file():
            raise FileNotFoundError("VCF file '%s' is not readable."%self.vcfpath)
        with open(self.vcfpath, 'r') as f:
            n = 0
            line = f.readline()
            #VCF header
            while line.startswith('##'):
                n += 1
                self.vcfheader.append(line)
                line = f.readline()
            if line.startswith(self.vcfHeaderExpected[0]):
                n += 1
                self.vcfheader.append(line)
            self.headerlinenumber = n
            #Check VCF formatting
            self.check_vcf_header(self.vcfheader[-1].split("\t"))
            
            #VCF body
            if not headerOnly:
                for line in f:
                    n += 1
                    self.vcfbody.append(line)
                self.totlines = n
    
    def count_lines(self):
        """Count the number of lines in a VCF file."""            
        if self.totlines == -1:
            if self.readable_file:
                self.totlines = sum(1 for _ in open(self.vcfpath))
            elif self.vcfbody != []:
                self.totlines = len(self.vcfbody)
            else:
                Config.error_print("Cannot assign number of lines")
        else:
            Config.verbose_print("Number of lines already calculated")
    
    def define_intervals(self, chunks):
        """
        Determine the evenly distributed ranges of lines in the VCF file matching each chunk

        chunks (int) : The number of chunks to divide the VCF file into
        
        return (list) : The unique sizes of all intervals

        """
        assert (self.totlines-self.headerlinenumber) >= chunks >= 1
        
        #To make the intervals as even as possible, they are formed using the formula:
        # n = i * (x+1) + j * x
        # with:
        # n = total number of variants,
        # x = n//(# of chunks),
        # i and j the numbers of intervals of size x or x+1 so that the sum equals n
        
        base_interval = (self.totlines - self.headerlinenumber)//chunks #corresponds to 'x'
        last_interval = self.totlines - self.headerlinenumber - (base_interval * (chunks - 1))
        extra = last_interval - base_interval #remainder is equal to 'i'
        majorshift = self.headerlinenumber #intervals must consider the VCF header as well

        self.intervals = [[majorshift+part*(base_interval+1),majorshift+(part+1)*(base_interval+1)] for part in range(extra)]
        self.intervals += [[majorshift+extra+part*(base_interval),majorshift+extra+(part+1)*(base_interval)] for part in range(extra,chunks)]
        
        uniquechunksizes = [size_inter for size_inter in set([interval[1]-interval[0] for interval in self.intervals])]
        if  min(uniquechunksizes) < 100:
            log = ("Low chunk size of %s detected,"
                    " consider increasing the chunk size upper limit or reducing the number of cores.")%min(uniquechunksizes)
            print(log)
        
        return uniquechunksizes
        
            

class Vcfchunk(Vcf):
    """Store variant chunks from a VCF file"""

    def __init__(self, chunknumber = 1, vcf = None):
        """
        chunknumber (int) : The rank of the chunk of variants to use
        vcf (Vcf) : the VCF object created from the corresponding VCF file
        vcfbody (list) : The list of lines in the VCF chunk

        """
        if vcf is not None:
            assert type(vcf).__name__ == "Vcf"
            self.__dict__.update(vcf.__dict__)
        
        if not self.readable_file():
            Config.error_print("The VCF file has to be read several times, it should not be spawned by a subprocess")
            raise FileNotFoundError("VCF file 'is not readable."%self.vcfpath)
        self.chunknumber = chunknumber
        self.vcfbody = []
    
    @property
    def interval(self):
        return [self.intervals[self.chunknumber-1][0],self.intervals[self.chunknumber-1][1]]
    
    def read_vcf(self):
        """Reads variants from a chunk of a VCF file"""
        firstline = max(self.headerlinenumber, self.interval[0])
        lastline = self.interval[1]
        self.vcfbody = []
        currentline = 0
        with open(self.vcfpath, 'r') as f:
            while currentline < lastline:
                currentline += 1
                if currentline > firstline:
                    self.vcfbody.append(f.readline())
                else:
                    f.readline()
