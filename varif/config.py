import argparse

class Config(object):
    options={}

    @staticmethod
    def set_options():
        parser = argparse.ArgumentParser(description='''
            Filter and annotate variants likely to be
            differentially expressed (or fixed) among your sample(s)
            wih varif.
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

        parser.add_argument('--fixed', dest='fixed', action='store_true',
        help='Add population-fixed variants in the result')
        parser.add_argument('--no-fixed', dest='fixed', action='store_false',
        help='Do not add population-fixed variants in the result')

        parser.add_argument('--all-variants', dest='allVariants', action='store_true',
        help='Add all variants in the result (including fixed)')
        parser.add_argument('--best-variants', dest='allVariants', action='store_false',
        help='Add only differentially expressed variants')

        parser.add_argument('--all-regions', dest='allRegions', action='store_true',
        help='Add variants from all regions of the genome')
        parser.add_argument('--gene-regions', dest='allRegions', action='store_false',
        help='Add variants only in gene-annotated regions')


        parser.add_argument('--nucl-window-before', dest='nuclWindowBefore', type=int, default=0,
        metavar="INT",
        help='Number of bases before the variant to be included')
        parser.add_argument('--nucl-window-after', dest='nuclWindowAfter', type=int, default=0,
        metavar="INT",
        help='Number of bases after the variant to be included')
        parser.add_argument('--prot-window-before', dest='protWindowBefore', type=int, default=0,
        metavar="INT",
        help='Number of amino acids before the variant to be included')
        parser.add_argument('--prot-window-after', dest='protWindowAfter', type=int, default=0,
        metavar="INT",
        help='Number of amino acids after the variant to be included')

        parser.add_argument('--depth', dest='mindepth', type=int, default=5,
        metavar="DEPTH",
        help='Minmal total read depth for a alt to be considered')
        parser.add_argument('--ratio-alt', dest='minaaf', type=float, default=0.8,
        metavar="AAF",
        help='Minmal ratio of alt/total depth to call it true alt')
        parser.add_argument('--ratio-no-alt', dest='maxraf', type=float, default=0.2,
        metavar="RAF",
        help='Maximal ratio of alt/total depth to call it true ref')
        parser.add_argument('--control', dest='refRanks', type=str, default="",
        metavar="IDi,...,IDn",
        help='Comma-separated list of negative control sample rank (order in VCF file, starting with 1)')

        parser.add_argument('--filtered-csv', dest='csv', type=str, default=None,
        metavar="FILE",
        help='Name of the output CSV file')
        parser.add_argument('--filtered-vcf', dest='filteredvcf', type=str, default=None,
        metavar="FILE",
        help='Name of the output VCF file')

        parser.add_argument('--version', dest='version', action='store_true',
        help='Show varif version')
        
        args=parser.parse_args()
        Config.options=vars(args)
