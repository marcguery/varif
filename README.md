# Varif
A filtering and annotating program for VCF-formatted variants.

# Installation and upgrade

```bash
git clone https://github.com/marcguery/varif.git
cd varif
pip3 install .
varif -h
varif --version
```

Or without using `pip`:

```bash
git clone https://github.com/marcguery/varif.git
export PYTHONPATH=$(pwd -P)/varif:$PYTHONPATH
./varif/bin/varif -h
./varif/bin/varif --version
```

# Inputs

All input files must have the same chromosome names.

## VCF (`-vcf`)

A VCF file formatted as described in `samtools` specifications at [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

Mandatory metadata headers are:

- #CHROM
- POS
- ID
- REF
- ALT
- QUAL
- FILTER
- INFO

## GFF (`-gff`)

A GFF3 file formatted as described in `ensembl` specifications at [https://www.ensembl.org/info/website/upload/gff3.html](https://www.ensembl.org/info/website/upload/gff3.html).

Supported types of features are:

- *gene* or containing the keyword 'gene' (*ncRNA_gene*, *protein_coding_gene*, *pseudogene*, ...)
  **Used to annotate the variants**.
- *CDS*
  **Used to translate sequences affected by the variants**.

## FASTA (`-fasta`)

The FASTA file containing one sequence per chromosome.

## PED (`--ped`)

The optional PED file containing the family membership and parental relationships between samples from the VCF file. This file should be formatted as described in https://zzz.bwh.harvard.edu/plink/data.shtml#ped.

Only samples present in the PED file will be processed.

# Overview

Use `varif` to filter and annotate variants submitted through a VCF file.

```
options:
  -h, --help                show this help message and exit
  -vcf FILE                 VCF file
  -gff FILE                 GFF3 file
  -fasta FILE               FASTA file
  -outfilename FILENAME     Name of the main output file (no extension) that
                            will be appended with the group names
  --ncores INT              Number of parallel jobs to run
  --chunk-size INT          Maximal number of variants to be processed in each
                            chunk
  --ped FILE                PED file
  --comparison STR          Compare variants between "families", "lineages",
                            "self" or "all"
  --nucl-window-before INT  Number of bases to include before the allele
  --nucl-window-after INT   Number of bases to include after the allele
  --prot-window-before INT  Number of amino acids to include before the allele
  --prot-window-after INT   Number of amino acids to include after the allele
  --depth INT               Minimal read depth for a sample variant to be
                            considered
  --ratio-alt FLOAT         Minimal sample proportion of allele/total read depth
                            to ignore other alleles
  --ratio-ref FLOAT         Maximal sample proportion of allele/total read depth
                            to ignore it
  --exclude-intergenic      Keep only variants from gene-annotated regions
  --max-heterozygous FLOAT  Maximal proportion of heterozygous loci in each
                            group
  --max-missing FLOAT       Maximal proportion of missing loci in each group
  --min-apf-diff FLOAT      Minimal difference in allele population frequency
                            between groups
  --min-samples-diff INT    Minimal difference in the number of samples with a
                            distinct allele between groups
  --min-maf1 FLOAT          Minimal population MAF reached in a group
  --max-maf1 FLOAT          Maximal population MAF reached in a group
  --min-maf2 FLOAT          Minimal population MAF reached in both groups
  --max-maf2 FLOAT          Maximal population MAF reached in both groups
  --output-vcf              Output filtered VCF file(s)
  --verbose                 Show more details
  --version                 Show varif version
```

For each allele in each sample, `varif` will calculate the allele sample proportion (see [Allele Sample Proportions](#asp)) and classify it as true homozygous or not based on filters applied with the available options. The values of all allele sample proportions are compared to find sites likely to be differentially mutated between group of samples.

Alleles passing the filters will be stored in a CSV file (multiple alleles at the same position are separated) generated with the name (without extension) set by the option `-outfilename`. If the option `--output-vcf` is set, the corresponding filtered variants (multiple alleles at the same position are merged) from the VCF file will be saved as well. Each combination of  groups compared will lead to a unique CSV (and VCF) file (see [Sample groups](#groups)).

Here is an example of the  `varif` main output obtained with default parameters:

| Chromosome | Position | Type                     | Ref  | Alt  | CDS_ref          | CDS_alt  | AA_ref     | AA_alt      | Annotation                   | Missingness_G1 | Missingness_G2 | Heterozygous_G1 | Heterozygous_G2 | Alt_frequency_G1 | Alt_frequency_G2 | S1       | S2       | S3       |
| ---------- | -------- | ------------------------ | ---- | ---- | ---------------- | -------- | ---------- | ----------- | ---------------------------- | -------------- | -------------- | --------------- | --------------- | ---------------- | ---------------- | -------- | -------- | -------- |
| Chr1       | 123123   | SNP                      | C    | G    | NA               | NA       | NA         | NA          | NA                           | 0.000000       | 0.000000       | 0.000000        | 0.000000        | 0.666667         | 0.666667         | 0.912643 | 0.000000 | 1.000000 |
| Chr1       | 456456   | SNP,gene,mRNA            | T    | A    | NA               | NA       | NA         | NA          | Unknown function (gene_id_1) | 0.000000       | 0.000000       | 0.333333        | 0.333333        | 0.333333         | 0.333333         | 0.540000 | 0.000000 | 1.000000 |
| Chr2       | 123123   | SNP,CDS,exon,gene,mRNA   | G    | T    | AAGCG (CDS_id_2) | AATCG    | L (45/458) | I (45/458)  | Unknown function (gene_id_2) | 0.333333       | 0.333333       | 0.500000        | 0.500000        | 0.500000         | 0.500000         | NA       | 0.000000 | 0.440000 |
| Chr3       | 123456   | INDEL,CDS,exon,gene,mRNA | G    | GATA | ATGAT (CDS_id_3) | ATGATAAT | D (52/123) | DN (52/124) | Unknown function (gene_id_3) | 0.000000       | 0.000000       | 0.000000        | 0.000000        | 0.000000         | 0.000000         | 0.000000 | 0.120000 | 0.000000 |

## CSV content

The output CSV file separated by semicolons contains:

- ***Chromosome*** : The chromosome name
- ***Position*** : The starting position of the allele
- ***Type*** : The type of allele (SNP, INDEL...) followed by genomic features found
- ***Ref*** : The reference sequence +/- up and downstream sequences (separated by '|')
- ***Alt*** : The alternate sequence
- ***CDS_ref*** : The windowed CDS reference sequence including the allele followed by the CDS identifier\*
- ***CDS_alt*** : The windowed CDS alternate sequences including the allele\*
- ***AA_ref*** : The amino acids resulting from the windowed reference sequence followed by the codon position\*
- ***AA_alt*** : The amino acids resulting from the windowed alternate sequences followed by the codon position\*
- ***Annotation*** : The annotation described in the GFF file followed by the gene identifier\*\*
- ***Missingness_G1*** : The proportion of group 1 samples with a missing locus
- ***Heterozygous_G1*** : The proportion of group 1 samples (excluding samples with missing locus) with a heterozygous loci
- ***Alt_frequency_G1*** : The alternate allele population frequency of group 1 samples (excluding samples with missing locus)
- The remaining columns contain the allele proportion of each sample (ASPs)

\*: Applies for CDS regions, otherwise *NA*

\*\*: Applies for transcribed regions, otherwise *NA*

Note that when the CDS_alt field is empty and that the _lt field contains '*(0/0)*', this means that both reference and mutated CDS sequences are identical (this can happen when an INDEL occurs at the edge of a CDS and an intron).

## <a name="asp"></a>Allele Sample Proportions (ASP)

There are 2 cut-offs of minimal alternate ASP (minVSP with `--ratio-alt`) and maximal reference ASP (maxRSP with `--ratio-ref`) to determine if a sample is carrying the reference or alternate allele. By default, minVSP is equal to 0.8 and maxRSP to 0.2.

For each allele and for each sample, the ASP is calculated using the different allele read depths (ARD):

***ASP*** = (*alternate ARD*) / (*all alternates ARD* + *reference ARD*)

If the ASP is above minVSP, the loci is considered homozygous for the alternate allele, while if the ASP is below maxRSP, the loci is considered homozygous for the non alternate allele. When the ASP is between maxRSP and minVSP, the genotype cannot be called and the variant is heterozygous.

## <a name="groups"></a>Sample groups

Samples can be separated into groups which can have their allele frequencies compared with the option `--comparison` with one of the keywords *families*,*self*, *lineages*, or *all* (the default). Samples are grouped into families or parental-offspring relationships according to the metadata from the PED file. Allele differential distribution can be compared:

- between every combination of families (*families*),
- within all members of each family (*self*),
- between direct parents and offspring  (*lineages*),
- or all-against-all VCF samples (*all*).

# Filters

## Minimal depth

The ASPs are calculated only if the sample read depth (the sum of all the reads REF and ALTs) is equal or superior to the value provided with the option `--depth`, with a default value of 5.

## Allele population frequencies

The options `--min-apf-diff` (default 0) and `--min-samples-diff` (default 0) compare the prevalence of homozygous alleles (determined for each sample from ASPs) between the groups defined by the PED file.

The option `--min-apf-diff` corresponds to the minimal difference in (either alternate or reference) allele frequencies between groups to keep a variant. Allele frequency is calculated with the formula:

***APF*** = (*homozygous samples for one allele*) / (*all homozygous samples for both alleles* + *heterozygous samples*)

The option `--min-samples-diff` refers to the number of samples with opposite alleles in each group. For example, a value of 2 will keep only variants with at least 2 samples with an alternate allele in a group and 2 samples with a reference allele in the other group. A value of 0 will not filter out any variant. The value of this option is limited by the number of samples in each group independently (potentially resulting in unequal cutoffs between groups).

## Genetic regions

By default, all genomic regions are included. However, alleles outside a gene (all of the reference bases are located outside of a gene or in an intron) can be removed from the output with the option `--exclude-intergenic`.

## Missing data

Variants with too many missing (because of insufficient depth) or heterozygous (ASP between maxRSP and minVSP) genotypes can be filtered out respectively with `--max-missing` and `--max-heterozygous` (default 1) which discards alleles if one group has more missing/heterozygous genotypes than the selected proportion.

## Population Minor Allele Frequencies (MAF)

Variants can be filtered out based on the value of the population MAF in one of the groups compared or both of them. 

A variant is kept if all these conditions are met:

- The MAFs of both groups are above `--min-maf2` (default 0)
- The MAF of one of the groups is above `--min-maf1` (default 0)
- The MAF of one of the groups is below `--max-maf1` (default 0.5)
- The MAFs of both groups are below `--max-maf2` (default 0.5)

Cutoffs of MAF should be set as follow: `--min-maf2` &le; `--min-maf1` &le; `--max-maf1` &le; `--max-maf2`.

# Sequence annotations

## Up and downstream genomic sequences

The genomic sequences that are upstream or downstream of the reference allele can be shown in the Ref column of the output CSV file with respectively `--nucl-window-before` and `--nucl-window-after`.

## Overlapping features

Alleles inside features are annotated given the feature type available in the third column of the GFF file. If the feature type contains the keyword 'gene', their annotation is saved in the Annotation column of the CSV file. If the feature type is 'CDS', the automatic translation will be processed and the CDS and amino acid sequences associated with the CDS annotation will be added in the CDS_ref, CDS_alt, AA_ref and AA_alt columns of the CSV file.
All the feature types detected are stored after the allele type in the Type column of the CSV file.

## Automatic translation

When the allele is inside a CDS, `varif` will predict the new CDS sequence with the alternate allele by replacing the reference sequence (without bases from intronic regions) by the mutation. However, knowing where to start the CDS can be difficult in the case of a mutation affecting the very first codon of the CDS. In that matter, several further modifications are added to produce the most biologically relevant CDS:

- When the mutated sequence is shorter than the reference sequence:
  - If the very first bases of the initial CDS are not recovered after adding the mutation, bases are added from the upstream sequence to lead to a mutated CDS of the same size than the initial one
  - Otherwise, bases from the downstream sequence are added to match the initial CDS size
- When the mutated sequence is longer than the reference sequence:
  - Bases located in the upstream region of the CDS are not included in the mutated CDS
  - All other bases from the mutation are kept (including those from the downstream region) and bases from the downstream sequence after the mutation are added to make the length of the mutated CDS a multiple of 3

Finally, reference bases located within introns are removed before being replaced by the mutated sequence. For now, as no base from mutated sequence are removed, it is possible that bases which are actually part of an intron are included in the mutated CDS.

Amino acids affected by the reference and alternate sequence are obtained by translating the initial and mutated CDS respectively. Additional codons can be added in a chosen window with `--prot-window-before` and `--prot-window-after` options; limited to the first and last codon (or when a stop codon is obtained). The position of the codon is then calculated from their rank from the beginning of the initial and mutated CDS and shown in the AA_ref and AA_alt columns of the CSV file. The codon positions are associated with the total number of codons in each CDS, including the stop codon.

# Optimization

## Multiprocessing and memory

The most time consuming step of the pipeline is when the ASPs are calculated. Variants are processed by chunk whose number is determined from an upper limit on the number of variants processed by chunk set by `--chunk-size` (1000 variants by chunk by default). The number of variants by chunk can be reduced (down to 100) to save up on memory, although this will increase the running time. Several chunks can be processed at the same time by increasing the number of cores with `--ncores` (1 by default). The size of chunks may be reduced so that the total number of chunks match the number of cores.

# Limitations

## Reference and alternate alleles

For now, `varif` compares allele frequencies of the alternate allele with allele frequencies of the **non-alternate** allele(s) rather than the reference allele. In this documentation, reference should then be read as non-alternate. Although confusing reference with non-alternate can be appropriate with bi-allelic variants, this is not the case with multi-allelic variants.

## Multi allelic CDS

As `varif` considers only one allele at a time, the wrong mutated CDS can be translated if two different mutations occur in the same one.

## CDS prediction accuracy

Since `varif` does not consider bases from a mutation if they reach the upstream sequence (before the first codon), a longer mutated CDS starting before the initial first codon could be missed.

Similarly, the end of the mutated CDS is not precisely determined but estimated by the addition of enough bases to reach at least the initial CDS size or a multiple of 3.

## Mutation overlapping CDS and introns

if an INDEL includes the edge between a CDS and an intron, the resulting protein could include false codons from the intronic regions as for now no mutated base within an intronic region are removed.

# Examples

1. 

   - Save all possible alleles regardless of their ASPs (calculated if &gt; 5 read depth)
   
   - Process in parallel 4 batches of 1000 variants
   
   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --output-vcf --chunk-size 1000 --ncores 4 \
       --depth 5 \
       --min-apf-diff 0 --min-samples-diff 0
   ```
   
2. Same as 1. +

   - Include only variants within genes
   
   - Add 10 bases upstream and downstream of the alleles
   
   - Add 5 amino acids upstream and downstream of the alleles
   
   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --output-vcf --chunk-size 1000 --ncores 4 \
       --depth 5 \
       --min-apf-diff 0 --min-samples-diff 0 \
       --exclude-intergenic \
       --nucl-window-before 10 --nucl-window-after 10 \
       --prot-window-before 5 --prot-window-after 5
   ```
   
3. Same as 2. +

   - Attribute alternate allele to samples with ASP &ge; 0.8, reference allele to samples with ASP &le; 0.2
   
   - Keep only variants differentially distributed among all samples : &ge; 1 sample with an allele and &ge; 1 sample with the other allele
   
   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --output-vcf --chunk-size 1000 --ncores 4 \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --min-apf-diff 0 --min-samples-diff 1 \
       --exclude-intergenic \
       --nucl-window-before 10 --nucl-window-after 10 \
       --prot-window-before 5 --prot-window-after 5 \
       --comparison all
   ```
   
4. Same as 2. +

   - Attribute alternate allele to samples with ASP &ge; 0.8, reference allele to samples with ASP &le; 0.2
   
   - Keep only variants differentially distributed between pairs of families of samples (PED file) : &ge; 1 sample with an allele in one family and &ge; 1 sample with the other allele in the other family
   
   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --output-vcf --chunk-size 1000 --ncores 4 \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --min-apf-diff 0 --min-samples-diff 1 \
       --exclude-intergenic \
       --nucl-window-before 10 --nucl-window-after 10 \
       --prot-window-before 5 --prot-window-after 5 \
       --ped my_ped.ped --comparison families
   ```
   
5. Same as 2. +

   - Attribute alternate allele to samples with ASP &ge; 0.8, reference allele to samples with ASP &le; 0.2
   
   - Keep only variants differentially distributed between pairs of parent/offspring of samples (PED file) : &ge; 2 samples with an allele in one parent/offspring and &ge; 2 samples with the other allele in the other parent/offspring
   
   - Keep variants with 80 % of samples with non-missing call in each group
   
   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --output-vcf --chunk-size 1000 --ncores 4 \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --min-apf-diff 0 --min-samples-diff 2 \
       --exclude-intergenic \
       --nucl-window-before 10 --nucl-window-after 10 \
       --prot-window-before 5 --prot-window-after 5 \
       --ped my_ped.ped --comparison lineages \
       --max-missing 0.2
   ```
