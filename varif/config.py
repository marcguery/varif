class Config(object):
    mandatoryOptions={
        "vcf":None,
        "gff":None,
        "fasta":None}
    intOptions={
        "mindepth":5,
        "maximum":99999, "minimum":-99999,
    }
    floatOptions={
        "maxprop":0.2, "minprop":0.8,
    }
    boolOptions={
        "fixed":False, "allVariants":False,
        "allRegions":False, "show":False,
    }
    strOptions={
        "csv":"Variants.csv"
    }
    options={**intOptions, **floatOptions, **boolOptions, **strOptions}
    lengthOptions=len(options)+len(mandatoryOptions)

    @staticmethod
    def set_options(options):
        for option in options:
            if option in Config.intOptions:
                try:
                    options[option]=int(options[option])
                except:
                    print("Integer option %s is not well formatted"%option)
                    continue
            elif option in Config.floatOptions:
                try:
                    options[option]=float(options[option])
                except:
                    print("Float option %s is not well formatted"%option)
                    continue
            elif option in Config.boolOptions:
                try:
                    assert options[option]==True
                except:
                    print("Integer option %s is not well formatted"%option)
                    continue
            elif option in Config.strOptions:
                try:
                    assert type(options[option]) is str
                except:
                    print("String option %s is not well formatted"%option)
                    continue               
            elif option in Config.mandatoryOptions:
                try:
                    assert type(options[option]) is str
                except:
                    print("String option %s is not well formatted"%option)
                    continue
            else:
                print("Ignored option %s"%option)
                continue
            
            Config.options[option]=options[option]
        if len(Config.options) < Config.lengthOptions:
            print("Some mandatory options:\n%s\nare missing from your options:\n%s"
            %(", ".join(Config.mandatoryOptions.keys()), ", ".join(options)))
            raise(Exception)