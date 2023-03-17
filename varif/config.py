import argparse

class Config(object):
    options={}

    @staticmethod
    def set_options():
        parser = argparse.ArgumentParser(description='''
            Filter and annotate alleles likely to be
            differentially mutated among samples by comparing 
            their Allele Sample Proportions (ASPs)
            ''', formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=30))
        
        parser.add_argument('-vcf', type=str, default=None,
        metavar="FILE",
        help='VCF file')

        parser.add_argument('-gff', type=str, default=None,
        metavar="FILE",
        help='GFF3 file')

        parser.add_argument('-fasta', type=str, default=None,
        metavar="FILE",
        help='FASTA file')
        
        parser.add_argument('-outfilename', dest='outFile', type=str, default=None,
        metavar="FILENAME",
        help='Name of the main output file (no extension) that will be appended with the group names')

        parser.add_argument('--ncores', dest='ncores', type=int, default=1,
        metavar="CORES",
        help='Number of parallel jobs to run')

        parser.add_argument('--ped', type=str, default=None,
        metavar="FILE",
        help='PED file')
        parser.add_argument('--comparison', dest='comparison', type=str, default="all",
        metavar="GROUP", help='Compare variants between "families", "lineages", "selfself" or "all"')

        parser.add_argument('--fixed', dest='fixed', action='store_true',
        help='Add population-fixed variants in the result')
        parser.add_argument('--no-fixed', dest='fixed', action='store_false',
        help='Do not add population-fixed mutations in the result')

        parser.add_argument('--all-variants', dest='allVariants', action='store_true',
        help='Add all variants in the result (including fixed mutations)')
        parser.add_argument('--best-variants', dest='allVariants', action='store_false',
        help='Add only differentially mutated alleles in the main output (the whole variant for the output VCF)')

        parser.add_argument('--all-regions', dest='allRegions', action='store_true',
        help='Add variants from all regions of the genome')
        parser.add_argument('--gene-regions', dest='allRegions', action='store_false',
        help='Add variants only in gene-annotated regions')


        parser.add_argument('--nucl-window-before', dest='nuclWindowBefore', type=int, default=0,
        metavar="INT",
        help='Number of bases to include before the allele')
        parser.add_argument('--nucl-window-after', dest='nuclWindowAfter', type=int, default=0,
        metavar="INT",
        help='Number of bases to include after the allele')
        parser.add_argument('--prot-window-before', dest='protWindowBefore', type=int, default=0,
        metavar="INT",
        help='Number of amino acids to include before the allele')
        parser.add_argument('--prot-window-after', dest='protWindowAfter', type=int, default=0,
        metavar="INT",
        help='Number of amino acids to include after the allele')

        parser.add_argument('--depth', dest='mindepth', type=int, default=5,
        metavar="DEPTH",
        help='Minimal read depth for a sample variant to be considered')
        parser.add_argument('--ratio-alt', dest='minaltasp', type=float, default=0.8,
        metavar="ASP",
        help='Minimal sample propotion of allele/total read depth to consider it fixed')
        parser.add_argument('--ratio-ref', dest='maxrefasp', type=float, default=0.2,
        metavar="ASP",
        help='Maximal sample propotion of allele/total read depth to ignore it')

        parser.add_argument('--max-missing', dest='maxMissing', type=float, default=1,
        metavar="PROP",
        help='Maximal proportion of missing or mixed ASP in each group')
        parser.add_argument('--max-similarity', dest='maxSimilarity', type=float, default=1,
        metavar="PROP",
        help='Maximal ratio of min(mutated samples prop)/max(mutated samples prop) between groups (non fixed mutations)')
        parser.add_argument('--min-mutated', dest='minMutations', type=float, default=0,
        metavar="PROP",
        help='Minimal proportion of mutated samples in the most mutated group (non fixed mutations)')

        parser.add_argument('--output-vcf', dest='outputVcf', action='store_true',
        help='Whether to output filtered VCF files or not')

        parser.add_argument('--version', dest='version', action='store_true',
        help='Show varif version')
        
        args=parser.parse_args()
        Config.options=vars(args)
