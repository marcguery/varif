import argparse
from sys import stderr
import os

class Config(object):
    """All variants found from a VCF file.
    
        options (dict): 
        long_options (dict): Descriptive names (values) of options (keys)
    """
    _options = {}
    options = {}
    long_options = {"vcf": "VCF file", "gff": "GFF file", "fasta": "FASTA file",
                               "outFile": "Output file name", "ncores": "Number of cores",
                               "chunksize": "Chunk size",
                               "ped": "PED file", "comparison": "Comparison",
                               "nuclWindowBefore": "Nucleotides before", "nuclWindowAfter": "Nucleotides after",
                               "protWindowBefore": "Aminoacids before", "protWindowAfter": "Aminoacids after",
                               "mindepth" : "Minimal depth", "minaltasp" : "Alt ASP", "maxrefasp": "Ref ASP", 
                               "intergenicRegions": "Intergenic regions", 
                               "maxHeteroz": "Max. ratio of heterozygous loci", "maxMissing": "Max. ratio of missing loci",
                               "minApfDiff": "Min. diff in APF", "minSamplesDiff": "Min. number of diff. samples",
                               "minMaf1": "Min. MAF (1 group)", "maxMaf1": "Max. MAF (1 group)",
                               "minMaf2": "Min. MAF (2 groups)", "maxMaf2": "Max. MAF (2 groups)",
                               "outputVcf": "Output VCF", "verbose":"Verbosity",
                               "version": "Version"}
    
    @staticmethod
    def error_print(*args, **kwargs):
        """Prints to STDERR"""
        print(*args, **kwargs, file = stderr)
        return
    
    @classmethod
    def verbose_print(cls, *args, **kwargs):
        """Prints to STDOUT only of the 'verbose' option was used"""
        if cls.options["verbose"]:
            print(*args, **kwargs)
        else:
            return
        
    @classmethod
    def check_options(cls):
        """
        Verify the integrity of the arguments passed from the command line

        """
        if cls._options["version"] is True:#If version asked, only save version in options
            cls.options = {"version":cls._options["version"]}
            return
        for arg in ["vcf", "gff", "fasta", "outFile"]:
            if cls._options[arg] is None:
                raise NameError("%s is required"%cls.long_options[arg])
        
        if cls._options["comparison"] is not None:
            allowed_values = ["families", "lineages", "self", "all"]
            if cls._options["comparison"] not in allowed_values:
                raise ValueError("%s can only be one of '%s'"%(cls.long_options["comparison"],"', '".join(allowed_values)))
            if cls._options["ped"] is None and cls._options["comparison"] != "all":
                raise NameError("%s is required when comparing groups ('%s' comparison)"%(cls.long_options["ped"],
                                                                                          cls._options["comparison"]))
        
        if cls._options["ncores"] < 1:
            raise ValueError("%s must be at least 1, not %s"%(cls.long_options["ncores"], cls._options["ncores"]))
        if cls._options["chunksize"] < 100:
            raise ValueError("%s must be at least 100, not %s"%(cls.long_options["chunksize"], cls._options["chunksize"]))
    
        for arg in ["nuclWindowBefore", "nuclWindowAfter", "protWindowBefore", "protWindowAfter"]:
            if cls._options[arg] < 0:
                raise ValueError("%s must be >= 0, not %s"%(cls.long_options[arg], cls._options[arg]))
        
        if cls._options["minSamplesDiff"] < 0:
            raise ValueError("%s must be >= 0, not %s"%(cls.long_options["minSamplesDiff"], cls._options["minSamplesDiff"]))
        
        for arg in ["minMaf1", "maxMaf1", "minMaf2", "maxMaf2"]:
            if not 0 <= cls._options[arg] <= 0.5:
                raise ValueError("%s shoud be between 0 and 0.5, not %s"%(cls.long_options[arg], cls._options[arg]))
        if cls._options["minMaf1"] < cls._options["minMaf2"] or cls._options["maxMaf1"] > cls._options["maxMaf2"]:
            cls._options["minMaf1"] = cls._options["minMaf2"] if cls._options["minMaf1"] < cls._options["minMaf2"] else cls._options["minMaf1"]
            cls._options["maxMaf1"] = cls._options["maxMaf2"] if cls._options["maxMaf1"] > cls._options["maxMaf2"] else cls._options["maxMaf1"]            
            print("Reassigned 'minMaf1' and/or 'maxMaf1' to custom values of 'minMaf2' and/or 'maxMaf2'")
        if not cls._options["minMaf2"] <= cls._options["minMaf1"] <= cls._options["maxMaf1"] <= cls._options["maxMaf2"]:
            raise ValueError("MAF cutoffs did not satisfy these conditions:"+
                             "\n %s (was %s) <= %s (was %s) <= %s (was %s) <= %s (was %s)"%(cls.long_options["minMaf2"],
                                                                                            cls._options["minMaf2"],
                                                                                                           cls.long_options["minMaf1"],
                                                                                                           cls._options["minMaf1"],
                                                                                                           cls.long_options["maxMaf1"],
                                                                                                           cls._options["maxMaf1"],
                                                                                                           cls.long_options["maxMaf2"],
                                                                                                           cls._options["maxMaf2"]))
        
        for arg in ["minaltasp", "maxrefasp", "maxMissing", "minApfDiff"]:
            if not 0 <= cls._options[arg] <= 1:
                raise ValueError("%s must be a proportion, was %s"%(cls.long_options[arg], cls._options[arg]))
        assert cls._options["maxrefasp"] < cls._options["minaltasp"], ("%s (was %s) should be <= to"
                                                                  " %s (was %s)")%(cls.long_options["maxrefasp"],
                                                                                   cls._options["maxrefasp"],
                                                                                   cls.long_options["minaltasp"],
                                                                                              cls._options["minaltasp"])
        
        if len(cls._options["outFile"].split("/")[-1]) > 100:
            raise ValueError("%s must not exceed 100 characters"%(cls.long_options["outFile"]))
        
        if not os.path.isdir(os.path.dirname(cls._options["outFile"])):
            raise FileNotFoundError("Directory '%s' does not exist"%(cls._options["outFile"]))
        
        cls.options = cls._options
    
    @classmethod
    def copy_options(cls, options):
        """Copy existing options

        Args:
            options (dict): Name of options (keys) and their value
        """
        cls._options = options
        cls.check_options()

    @classmethod
    def load_options(cls):
        """Reads options from the command line using argparse"""
        parser = argparse.ArgumentParser(description = '''
            Filter and annotate alleles likely to be
            differentially mutated among samples by comparing 
            their Allele Sample Proportions (ASPs)
            ''', formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 30, width = 80))
        
        parser.add_argument('-vcf', type = str, default = None,
                            metavar = "FILE",
                            help = 'VCF file')

        parser.add_argument('-gff', type = str, default = None,
                            metavar = "FILE",
                            help = 'GFF3 file')

        parser.add_argument('-fasta', type = str, default = None,
                            metavar = "FILE",
                            help = 'FASTA file')
        
        parser.add_argument('-outfilename', dest = 'outFile', type = str, default = None,
                            metavar = "FILENAME",
                            help = 'Name of the main output file (no extension) that will be appended with the group names')

        parser.add_argument('--ncores', dest = 'ncores', type = int, default = 1,
                            metavar = "INT",
                            help = 'Number of parallel jobs to run')
        parser.add_argument('--chunk-size', dest = 'chunksize', type = int, default = 1000,
                            metavar = "INT",
                            help = 'Maximal number of variants to be processed in each chunk')

        parser.add_argument('--ped', type = str, default = None,
                            metavar = "FILE",
                            help = 'PED file')
        parser.add_argument('--comparison', dest = 'comparison', type = str, default = "all",
                            metavar = "STR", 
                            help = 'Compare variants between "families", "lineages", "self" or "all"')

        parser.add_argument('--nucl-window-before', dest = 'nuclWindowBefore', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of bases to include before the allele')
        parser.add_argument('--nucl-window-after', dest = 'nuclWindowAfter', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of bases to include after the allele')
        parser.add_argument('--prot-window-before', dest = 'protWindowBefore', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of amino acids to include before the allele')
        parser.add_argument('--prot-window-after', dest = 'protWindowAfter', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of amino acids to include after the allele')

        parser.add_argument('--depth', dest = 'mindepth', type = int, default = 5,
                            metavar = "INT",
                            help = 'Minimal read depth for a sample variant to be considered')
        parser.add_argument('--ratio-alt', dest = 'minaltasp', type = float, default = 0.8,
                            metavar = "FLOAT",
                            help = 'Minimal sample proportion of allele/total read depth to ignore other alleles')
        parser.add_argument('--ratio-ref', dest = 'maxrefasp', type = float, default = 0.2,
                            metavar = "FLOAT",
                            help = 'Maximal sample proportion of allele/total read depth to ignore it')
        
        parser.add_argument('--exclude-intergenic', dest = 'intergenicRegions', action = 'store_false',
                            default = True,
                            help = 'Keep only variants from gene-annotated regions')
        parser.add_argument('--max-heterozygous', dest = 'maxHeteroz', type = float, default = 1,
                            metavar = "FLOAT",
                            help = 'Maximal proportion of heterozygous loci in each group')
        parser.add_argument('--max-missing', dest = 'maxMissing', type = float, default = 1,
                            metavar = "FLOAT",
                            help = 'Maximal proportion of missing loci in each group')
        parser.add_argument('--min-apf-diff', dest = 'minApfDiff', type = float, default = 0,
                            metavar = "FLOAT",
                            help = 'Minimal difference in allele population frequency between groups')
        parser.add_argument('--min-samples-diff', dest = 'minSamplesDiff', type = int, default = 0,
                            metavar = "INT",
                            help = 'Minimal difference in the number of samples with a distinct allele between groups')
        
        parser.add_argument('--min-maf1', dest = 'minMaf1', type = float, default = 0,
                            metavar = "FLOAT",
                            help = 'Minimal population MAF reached in a group')
        parser.add_argument('--max-maf1', dest = 'maxMaf1', type = float, default = 0.5,
                            metavar = "FLOAT",
                            help = 'Maximal population MAF reached in a group')
        parser.add_argument('--min-maf2', dest = 'minMaf2', type = float, default = 0,
                            metavar = "FLOAT",
                            help = 'Minimal population MAF reached in both groups')
        parser.add_argument('--max-maf2', dest = 'maxMaf2', type = float, default = 0.5,
                            metavar = "FLOAT",
                            help = 'Maximal population MAF reached in both groups')

        parser.add_argument('--output-vcf', dest = 'outputVcf', action = 'store_true',
                            help = 'Output filtered VCF file(s)')
        
        parser.add_argument('--verbose', dest = 'verbose', action = 'store_true',
                            help = 'Show more details')

        parser.add_argument('--version', dest = 'version', action = 'store_true',
                            help = 'Show varif version')
        
        args = parser.parse_args()
        cls._options = vars(args)
        cls.check_options()
            
