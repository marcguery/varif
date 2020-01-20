#!/usr/bin/env python3
import math
import vcfpy

class Varif(object):

    def frac_AD(recordSample):
        """
        Gets the ratios of each allele with ref allele.
        Note that 0/0 is 'considered' equal to 1.

        recordSample (list) : The ref allele followed by the alternate allele(s)

        return (list) : the ratios
        """
        #Handling division by zero, when there is no ref
        if sum(recordSample)==0:
            recordSample=[math.nan]*len(recordSample)
        else:
            recordSample=list(map(lambda x : x/sum(recordSample), recordSample))
        return recordSample

    def read_gff(gff, annotation):
        """
        Reads a GFF file

        gff (str) : the gff path location
        annotation (str) : annotation to be taken to subset positions

        return (dict) : The positions (start and end), ids and functions
        of extracted information per chromosome
        """
        gffInfo={}
        with open(gff) as f:
            rank=0
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                lineSplitted=line.split("\t", 8)
                if lineSplitted[2]==annotation:
                    rank+=1
                    identifier=lineSplitted[-1].split(";")[0][3:]
                    description=lineSplitted[-1].split(";")[-1][12:].rstrip("\n") if lineSplitted[-1].split(";")[-1].startswith('description') else 'Unknown function'
                    rangePosition=[int(lineSplitted[3]), int(lineSplitted[4])+1]
                    if not lineSplitted[0] in gffInfo:
                        gffInfo[lineSplitted[0]]={}
                        gffInfo[lineSplitted[0]]['positions']=[(rank, rangePosition)]
                        gffInfo[lineSplitted[0]]['ids']=[(rank, identifier)]
                        gffInfo[lineSplitted[0]]['descriptions']=[(rank, description)]
                    else:
                        gffInfo[lineSplitted[0]]['positions'].append((rank, rangePosition))
                        gffInfo[lineSplitted[0]]['ids'].append((rank, identifier))
                        gffInfo[lineSplitted[0]]['descriptions'].append((rank, description))
        for chromosome in gffInfo:
            gffInfo[chromosome]['positions']=sorted(gffInfo[chromosome]['positions'], key=lambda x : x[1][0])
        return gffInfo

    def get_id_from_position(position, gffreader, chromosome):
        """
        Retrieve information about an annotation level in which there is a variant

        position (int) : the position of the variant
        gffreader (read_gff) : the gff reader object
        chromosome (str) : the chromosome name in gff reader

        return (list) : ID, function, start/end positions or None if outside of a variant
        """
        finished=False
        maxValue=len(gffreader[chromosome]['positions'])-1
        minValue=0
        index=math.floor(maxValue-minValue/2)
        while not finished:
            mini=gffreader[chromosome]['positions'][index][1][0]
            maxi=gffreader[chromosome]['positions'][index][1][1]
            if position >= mini and position <= maxi:
                finished=True
                rank=0
                while gffreader[chromosome]['ids'][rank][0] != gffreader[chromosome]['positions'][index][0]:
                    rank+=1
                identifier=gffreader[chromosome]['ids'][rank][1]
                description=gffreader[chromosome]['descriptions'][rank][1]
                return [identifier, description, gffreader[chromosome]['positions'][index]]
            elif position > maxi and index < maxValue:
                minValue=index+1
                index=math.ceil((maxValue+index)/2)
            elif position < mini and index > minValue:
                maxValue=index-1
                index=math.floor((index+minValue)/2)
            else:
                finished=True

        return ["No ID", "No function", "No Positions"]

    def read_vcf(vcf):
        """Reads a VCF file with vcfpy library"""
        return vcfpy.Reader.from_path(vcf)

    def get_ratios_from_vcf(vcfreader, samples, mindepth=5):
        """
        Retrieve the allele ratios of variants in a vcf reader object

        vcfreader (VCF reader) : The VCF reader object
        mindepth : The minimal total depth to consider a variant
        """
        ratios={}
        for record in vcfreader:
            recordname=record.CHROM+":"+str(record.POS)
            ref=record.REF
            alt=record.ALT
            id=1
            while recordname in ratios:
                recordname+="."+str(id)
                id+=1
            ratios[recordname]={}
            recordsAD={call.sample:call.data.get('AD') for call in record.calls}
            for key in samples:
                ad=recordsAD[samples[key]] if sum(recordsAD[samples[key]])>mindepth else [0]*len(recordsAD[samples[key]])
                ratio=Varif.frac_AD(ad)
                ratios[recordname][key]=ratio
        return ratios

    def get_score_from_ratios(ratios, maxprop=0.1, minprop=None):
        minprop=1-maxprop if minprop is None else minprop
        assert 0 < maxprop < 0.5 and 0.5 < minprop < 1
        varState={sample:ratios[sample] > minprop for sample in ratios}
        minratio=min(ratios.values())
        maxratio=max(ratios.values())
        if math.isnan(minratio) or math.isnan(maxratio):
            score=math.nan
        #fixedVariants
        elif minratio > minprop:
            score = 0
        #Ambiguous variants
        elif minratio >= maxprop or maxratio <= minprop:
            score = -maxratio/minratio if minratio != 0 else -math.inf
        #True Variants
        elif minratio < maxprop and maxratio > minprop:
            score = maxratio/minratio if minratio != 0 else math.inf
        return [score, varState]


    def write_scores(scores, variants, gff, annotation, csv="Variants.csv", show=False, fixedvariants=True, allregions=False, allvariants=False):
        cutoff=1 if allvariants is False else -math.inf
        sortedKeys=sorted(scores, key= lambda x : (x.split(":")[0], int(x.split(":")[1].split(".")[0])))
        sortedvariantsKeys=sorted(variants[next(iter(variants))])
        with open(csv, 'w') as variantFile:
            variantFile.write(",".join(["Location", annotation, "Function", "Score"]+sortedvariantsKeys)+"\n")
            for key in sortedKeys:
                if any(score == 0 and fixedvariants is True or score >= cutoff for score in scores[key]):
                    identifier=Varif.get_id_from_position(int(key.split(":")[1].split(".")[0]), gff, key.split(":")[0])
                    if identifier[0] != "No ID" or allregions is True:
                        line=key+","+str(identifier[0])+","+str(identifier[1])+","+":".join([str(score) for score in scores[key]])+","+",".join([":".join(map(str, variants[key][sample])) for sample in sortedvariantsKeys])
                        if show is True:
                            print(line)
                        variantFile.write(line+"\n")

    
