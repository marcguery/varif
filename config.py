class Config(object):
    options={
        "mindepth":5, 
        "maxprop":0.2, "minprop":0.8,
        "maximum":99999, "minimum":-99999,
        "fixed":True, "allVariants":False,
        "allRegions":False, "show":False,
        "csv":"Variants.csv"

    }
    mandatoryOptions={
        "vcf":None,
        "gff":None,
        "fasta":None}
    lengthOptions=len(options)+len(mandatoryOptions)
    @staticmethod
    def set_options(options):
        for option in options:
            if option not in Config.options and option not in Config.mandatoryOptions:
                print("Ignored option %s"%option)
                continue
            Config.options[option]=options[option]
        if len(Config.options) < Config.lengthOptions:
            print("Some mandatory options:\n%s\nare missing from your options:\n%s"
            %(", ".join(Config.mandatoryOptions.keys()), ", ".join(options)))
            raise(Exception)