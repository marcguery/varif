import argparse
from sys import stderr
import os
from .version import __version__

class Config(object):
    """Varif common functions and options.
    
        options (dict): Options (keys) and their value
        long_options (dict): Descriptive names (values) of options (keys)
    """
    version = __version__
    _options = {}
    options = {}
    long_options = {"vcf": "VCF file", "gff": "GFF file", "fasta": "FASTA file",
                               "outFile": "Output file name", "ncores": "Number of cores",
                               "chunksize": "Chunk size",
                               "ped": "PED file", "comparison": "Comparison",
                               "nuclWindowBefore": "Nucleotides before", "nuclWindowAfter": "Nucleotides after",
                               "protWindowBefore": "Aminoacids before", "protWindowAfter": "Aminoacids after",
                               "mindepth" : "Minimal depth", 
                               "minVafHomozygous" : "Min. VAF for homozygous variants", "minVafPresent": "Min. VAF for present variants", 
                               "excludeIntergenic": "Exclude intergenic regions", 
                               "maxHeteroz": "Max. ratio of heterozygous loci", "maxMissing": "Max. ratio of missing loci",
                               "minApfDiff": "Min. diff in APF", "minSamplesDiff": "Min. number of diff. samples",
                               "minMaf1": "Min. MAF (1 group)", "maxMaf1": "Max. MAF (1 group)",
                               "minMaf2": "Min. MAF (2 groups)", "maxMaf2": "Max. MAF (2 groups)",
                               "outputVcf": "Output VCF", "verbose":"Verbose",
                               "version": "Version"}
    
    @staticmethod
    def error_print(*args, **kwargs):
        """Prints to STDERR"""
        print("[varif error]", *args, **kwargs, file = stderr)
        return
    
    @classmethod
    def verbose_print(cls, *args, **kwargs):
        """Prints to STDOUT only of the 'verbose' option was used"""            
        if cls.options["verbose"]:
            for arg in args:
                for splitline in arg.rstrip("\n").split("\n"):
                    print("[varif info]", splitline, **kwargs)
        else:
            return
    
    @classmethod
    def print_options(cls):
        """Prints all the options and their attributed value
        """
        longest_name = max(len(cls.long_options[opt]) for opt in cls.long_options)
        print("Running Varif %s with options:"%cls.version)
        for opt in cls.long_options:
            print("    {name:>{width}s} | {value}".format(name = cls.long_options[opt],
                                                     width = longest_name,
                                                     value = cls.options[opt]))
        
    @classmethod
    def check_options(cls):
        """
        Verify the integrity of the arguments passed from the command line

        """
        holdprint = []
        for arg in cls.long_options:
            if arg not in cls._options:
                cls.error_print("Missing documented option '%s' in command-line"%arg)
                raise KeyError
            
        if sorted(cls.long_options.keys()) != sorted(cls._options.keys()):
            suppoptions = set(cls._options.keys()).difference(cls.long_options.keys())
            cls.error_print("Options '%s' are not documented"%("', '".join(suppoptions)))
            raise KeyError()
                
        if cls._options["version"] is True:#If version asked, only save version in options
            cls.options = {"version":cls._options["version"]}
            return
        for arg in ["vcf", "gff", "fasta", "outFile"]:
            if cls._options[arg] is None:
                cls.error_print("%s is required"%cls.long_options[arg])
                raise NameError()
        
        if cls._options["comparison"] is not None:
            allowed_values = ["families", "lineages", "self", "all"]
            if cls._options["comparison"] not in allowed_values:
                cls.error_print("%s can only be one of '%s'"%(cls.long_options["comparison"],"', '".join(allowed_values)))
                raise ValueError()
            if cls._options["ped"] is None and cls._options["comparison"] != "all":
                cls.error_print("%s is required when comparing groups ('%s' comparison)"%(cls.long_options["ped"],
                                                                                          cls._options["comparison"]))
                raise NameError()
        
        if cls._options["ncores"] < 1:
            cls.error_print("%s must be at least 1, not %s"%(cls.long_options["ncores"], cls._options["ncores"]))
            raise ValueError()
        if cls._options["chunksize"] < 100:
            cls.error_print("%s must be at least 100, not %s"%(cls.long_options["chunksize"], cls._options["chunksize"]))
            raise ValueError()
    
        for arg in ["nuclWindowBefore", "nuclWindowAfter", "protWindowBefore", "protWindowAfter"]:
            if cls._options[arg] < 0:
                cls.error_print("%s must be >= 0, not %s"%(cls.long_options[arg], cls._options[arg]))
                raise ValueError()
        
        if cls._options["minSamplesDiff"] < 0:
            cls.error_print("%s must be >= 0, not %s"%(cls.long_options["minSamplesDiff"], cls._options["minSamplesDiff"]))
            raise ValueError()
        
        for arg in ["minMaf1", "maxMaf1", "minMaf2", "maxMaf2"]:
            if not 0 <= cls._options[arg] <= 0.5:
                cls.error_print("%s shoud be between 0 and 0.5, not %s"%(cls.long_options[arg], cls._options[arg]))
                raise ValueError()
        if cls._options["minMaf1"] < cls._options["minMaf2"]:
            oldminmaf1 = cls._options["minMaf1"]
            cls._options["minMaf1"] = cls._options["minMaf2"]
            holdprint.append("Reassigned '%s' (was %s) to the value of '%s' (%s)"%(cls.long_options["minMaf1"],
                                                                                   oldminmaf1,
                                                                                   cls.long_options["minMaf2"],
                                                                                   cls._options["minMaf2"]))
        if cls._options["maxMaf1"] > cls._options["maxMaf2"]:
            oldmaxmaf1 = cls._options["maxMaf1"]
            cls._options["maxMaf1"] = cls._options["maxMaf2"]
            holdprint.append("Reassigned '%s' (was %s) to the value of '%s' (%s)"%(cls.long_options["maxMaf1"],
                                                                                   oldmaxmaf1,
                                                                                   cls.long_options["maxMaf2"],
                                                                                   cls._options["maxMaf2"]))
        if not cls._options["minMaf2"] <= cls._options["minMaf1"] <= cls._options["maxMaf1"] <= cls._options["maxMaf2"]:
            cls.error_print("MAF cutoffs did not satisfy these conditions:"+
                             "\n %s (was %s) <= %s (was %s) <= %s (was %s) <= %s (was %s)"%(cls.long_options["minMaf2"],
                                                                                            cls._options["minMaf2"],
                                                                                            cls.long_options["minMaf1"],
                                                                                            cls._options["minMaf1"],
                                                                                            cls.long_options["maxMaf1"],
                                                                                            cls._options["maxMaf1"],
                                                                                            cls.long_options["maxMaf2"],
                                                                                            cls._options["maxMaf2"]))
            raise ValueError()
        
        if cls._options["minVafPresent"] is None:
            cls._options["minVafPresent"] = cls._options["minVafHomozygous"]
            holdprint.append("Assigned %s to the value of %s (%s)"%(cls.long_options["minVafPresent"], 
                                                                    cls.long_options["minVafHomozygous"],
                                                                    cls._options["minVafHomozygous"]))
        
        for arg in ["minVafHomozygous", "minVafPresent", "maxMissing", "minApfDiff"]:
            if not 0 <= cls._options[arg] <= 1:
                cls.error_print("%s must be a proportion, was %s"%(cls.long_options[arg], cls._options[arg]))
                raise ValueError()
        if cls._options["minVafPresent"] > cls._options["minVafHomozygous"]:
            cls.error_print("%s (was %s) should be <= to %s (was %s)"%(cls.long_options["minVafPresent"],
                                                                       cls._options["minVafPresent"],
                                                                       cls.long_options["minVafHomozygous"],
                                                                       cls._options["minVafHomozygous"]))
            raise ValueError()
        
        if len(cls._options["outFile"].split("/")[-1]) > 100:
            cls.error_print("%s must not exceed 100 characters"%(cls.long_options["outFile"]))
            raise ValueError()
        
        if os.path.dirname(cls._options["outFile"]) != "" and not os.path.isdir(os.path.dirname(cls._options["outFile"])):
            cls.error_print("Directory '%s' does not exist"%(cls._options["outFile"]))
            raise FileNotFoundError()
        
        cls.options = cls._options
        for material in holdprint:
            cls.verbose_print(material)
    
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
            their Variant Allele Frequencies (VAFs)
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
                            help = 'Number of parallel jobs to run [%(default)i]')
        parser.add_argument('--chunk-size', dest = 'chunksize', type = int, default = 1000,
                            metavar = "INT",
                            help = 'Maximal number of variants to be processed in each chunk [%(default)i]')

        parser.add_argument('--ped', type = str, default = None,
                            metavar = "FILE",
                            help = 'PED file')
        parser.add_argument('--comparison', dest = 'comparison', type = str, default = "all",
                            metavar = "STR", 
                            help = 'Compare variants between "families", "lineages", "self" or "all" [%(default)s]')

        parser.add_argument('--nucl-window-before', dest = 'nuclWindowBefore', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of bases to include before the allele [%(default)i]')
        parser.add_argument('--nucl-window-after', dest = 'nuclWindowAfter', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of bases to include after the allele [%(default)i]')
        parser.add_argument('--prot-window-before', dest = 'protWindowBefore', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of amino acids to include before the allele [%(default)i]')
        parser.add_argument('--prot-window-after', dest = 'protWindowAfter', type = int, default = 0,
                            metavar = "INT",
                            help = 'Number of amino acids to include after the allele [%(default)i]')

        parser.add_argument('--depth', dest = 'mindepth', type = int, default = 5,
                            metavar = "INT",
                            help = 'Minimal read depth for a sample variant to be considered [%(default)i]')
        parser.add_argument('--min-vaf-homozygous', dest = 'minVafHomozygous', type = float, default = 0.8,
                            metavar = "FLOAT",
                            help = 'Minimal VAF for a variant to be considered homozygous [%(default).2f]')
        parser.add_argument('--min-vaf-present', dest = 'minVafPresent', type = float, default = 0.05,
                            metavar = "FLOAT",
                            help = 'Minimal VAF for a variant to be considered present [%(default).2f]')
        
        parser.add_argument('--exclude-intergenic', dest = 'excludeIntergenic', action = 'store_true',
                            default = False,
                            help = 'Exclude variants outside of gene-annotated regions')
        parser.add_argument('--max-heterozygous', dest = 'maxHeteroz', type = float, default = 1,
                            metavar = "FLOAT",
                            help = 'Maximal proportion of heterozygous loci in each group [%(default).2f]')
        parser.add_argument('--max-missing', dest = 'maxMissing', type = float, default = 1,
                            metavar = "FLOAT",
                            help = 'Maximal proportion of missing loci in each group [%(default).2f]')
        parser.add_argument('--min-apf-diff', dest = 'minApfDiff', type = float, default = 0,
                            metavar = "FLOAT",
                            help = 'Minimal difference in allele population frequency between groups [%(default).2f]')
        parser.add_argument('--min-samples-diff', dest = 'minSamplesDiff', type = int, default = 0,
                            metavar = "INT",
                            help = 'Minimal difference in the number of samples with a distinct allele between groups [%(default)i]')
        
        parser.add_argument('--min-maf1', dest = 'minMaf1', type = float, default = 0,
                            metavar = "FLOAT",
                            help = 'Minimal population MAF reached in a group [%(default).2f]')
        parser.add_argument('--max-maf1', dest = 'maxMaf1', type = float, default = 0.5,
                            metavar = "FLOAT",
                            help = 'Maximal population MAF reached in a group [%(default).2f]')
        parser.add_argument('--min-maf2', dest = 'minMaf2', type = float, default = 0,
                            metavar = "FLOAT",
                            help = 'Minimal population MAF reached in both groups [%(default).2f]')
        parser.add_argument('--max-maf2', dest = 'maxMaf2', type = float, default = 0.5,
                            metavar = "FLOAT",
                            help = 'Maximal population MAF reached in both groups [%(default).2f]')

        parser.add_argument('--output-vcf', dest = 'outputVcf', action = 'store_true',
                            help = 'Output filtered VCF file(s) [%(default)s]')
        
        parser.add_argument('--verbose', dest = 'verbose', action = 'store_true',
                            help = 'Show more details')

        parser.add_argument('--version', dest = 'version', action = 'store_true',
                            help = 'Show varif version')
        
        args = parser.parse_args()
        cls._options = vars(args)
        cls.check_options()
            
