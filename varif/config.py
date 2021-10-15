import argparse

class Config(object):
    options={}

    @staticmethod
    def set_options():
        parser=argparse.ArgumentParser(
            description='''
            Filter and annotate variants likely to be
            differentially expressed (or fixed) among your sample(s)
            wih varif.
            ''')
        parser.add_argument('-vcf', required=True, type=str,
        metavar="FILE",
        help='VCF file')

        parser.add_argument('-gff', required=True, type=str,
        metavar="FILE",
        help='GFF3 file')

        parser.add_argument('-fasta', required=True, type=str,
        metavar="FILE",
        help='FASTA file')

        parser.add_argument('--fixed', dest='fixed', action='store_true',
        help='Add fixed variants in the result')
        parser.add_argument('--no-fixed', dest='fixed', action='store_false',
        help='Do not add fixed variants in the result')

        parser.add_argument('--all-variants', dest='allVariants', action='store_true',
        help='Add all variants in the result')
        parser.add_argument('--best-variants', dest='allVariants', action='store_false',
        help='Add only variants whose score is positive')

        parser.add_argument('--all-regions', dest='allRegions', action='store_true',
        help='Add variants from all regions of the genome')
        parser.add_argument('--gene-regions', dest='allRegions', action='store_false',
        help='Add variants only in gene-annotated regions')


        parser.add_argument('--window-before', dest='windowBefore', type=int, default=2,
        metavar="INT",
        help='Number of bases before the variant to be included')
        parser.add_argument('--window-after', dest='windowAfter', type=int, default=2,
        metavar="INT",
        help='Number of bases before the variant to be included')
        

        parser.add_argument('--show', dest='show', action='store_true',
        help='Print filtered VCF to stdout')
        parser.add_argument('--no-show', dest='show', action='store_false',
        help='Do not print any result')

        parser.add_argument('--depth', dest='mindepth', type=int, default=5,
        metavar="DEPTH",
        help='Minmal total read depth for a alt to be considered')
        parser.add_argument('--ratio-alt', dest='minprop', type=float, default=0.8,
        metavar="RATIO",
        help='Minmal ratio of alt/total depth to call it true alt')
        parser.add_argument('--ratio-no-alt', dest='maxprop', type=float, default=0.2,
        metavar="RATIO",
        help='Maximal ratio of alt/total depth to call it true ref')

        parser.add_argument('--csv', dest='csv', type=str, default=None,
        metavar="FILE",
        help='Name of the CSV file to be written')
        parser.add_argument('--filteredvcf', dest='filteredvcf', type=str, default=None,
        metavar="FILE",
        help='Name of the VCF file to be written')

        parser.add_argument('--max-score', dest='maximum', type=int, default=99999,
        metavar="SCORE",
        help='Maximal score attributed to positive-score variants')
        parser.add_argument('--min-score', dest='minimum', type=int, default=-99999,
        metavar="SCORE",
        help='Minimal score attributed to negative-score variants')
        args=parser.parse_args()
        Config.options=vars(args)
