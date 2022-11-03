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

# Summary

Use `varif` to filter and rank variants submitted through a VCF file.

For each variant in each sample, `varif` will calculate the Variant Alternate Frequency (VAF) and classify it as a true variant or not based on filters applied with the available options. The values of VAFs for each variant are compared to find sites likely to be differentially distributed among samples.

# Output

Variants passing the filters will be stored in a CSV file (variants from the same chromosome location are separated) generated with the name (without extension) set by the option `--out-filename`. If the option `--output-vcf` is set, the corresponding filtered variants from the VCF file will be saved as well. Each combination of  groups compared will lead to a unique CSV (and VCF) file (see [Sample groups](#groups)).

## CSV content

The output CSV file separated by semicolons contains:

- Chromosome: The chromosome name
- Position: The starting position of the variant
- Type: The type of variant (SNP, INDEL...) followed by its genomic location
- Ref: The reference sequence +/- up and downstream sequences (separated by '|')
- Alt: The alternate sequence
- CDSref: The windowed CDS reference sequence including the variant followed by the CDS identifier\*
- CDSalt: The windowed CDS alternate sequences including the variant\*
- AAref: The amino acids resulting from the windowed reference sequence followed by the codon position\*
- AAalt: The amino acids resulting from the windowed alternate sequences followed by the codon position\*
- Annotation: The annotation described in the GFF file followed by the gene identifier\*\*
- Score: The score of the variant (see [Variant scores](#scores))
- Sample#1: Proportion of alternate allele for sample 1
- ...
- Sample#n: Proportion of alternate allele for sample n

\*: Applies for CDS regions, otherwise *NA*

\*\*: Applies for transcribed regions, otherwise *NA*

# Filters

## Minimal depth

The variant allele frequencies are calculated only if the sample allele depth (the sum of all the reads REF and ALTs) is equal or superior to the value provided with the option `--depth`, with a default value of 5.

## Allele frequencies

There are 2 cut-offs of minimal Alternate Allele Frequency (minAAF with `--ratio-alt`) and maximal Reference Allele Frequency (maxRAF with `--ratio-no-alt`) used by `varif` to determine if variants are differentially present in a population. By default, minAAF is equal to 0.8 and maxRAF to 0.2.

For each variant at a chromosomal location and for each sample, the variant allele frequency (VAF) is calculated using the different allele depths (AD):

***VAF*** = (*variant AD*) / (*other variants AD* + *reference AD*)

If the VAF is above minAAF, the variant is considered a true variant, while if the VAF is below maxRAF, it is considered a true reference. When the VAF is between maxRAF and minAAF, the genotype cannot be called and the variant is mixed.

## <a name="groups"></a>Sample groups

By default, samples are not grouped into subcategories which implies that a self-self comparison will be made to determine which variant is differentially prevalent based on the presence of true variant samples and true reference samples. However, it is possible to separate samples into groups and look for variants differentially present between those groups by providing the option `--comparison` with one of the keywords *families*, *lineages* or *both*. Samples can be grouped either into families, or parental-offspring relationships according to the metadata from the PED file. Variant differential prevalence will be compared in every combination of families and direct parent versus offspring. The relationship between the two groups is symmetrical, meaning that variants detected in either group are candidates for being considered differentially present.

Groups with too many missing genotypes (because of insufficient depth or a mixed variant) can be filtered out with `--max-missing` (default 1) which discards variants if one group has more missing genotypes than the selected proportion.

Group-specific variants can be selected with the option `--max-similarity` (default 1) by providing the maximal proportion of the number of mutated samples between the group with less mutated samples over the one with the most mutated samples.

Variants that are too rare can be filtered out by the option `--min-variants` (default 0) which will discard variants that are present in less than this proportion of the number of samples in both group.

The combination of VAFs are used to classify variants that are:

- Differentially distributed between two groups:
  At least one VAF of one group is a true reference (below or equal to maxRAF) and at least one VAF of the other group is a true variant (above or equal to minAAF).
- Fixed in the population:
  All VAFs of the samples are above or equal to minAAF or below or equal to maxRAF. Either the option `--all-variants` or `--fixed` will show these variants.
- Not differentially distributed among samples, for all other cases. The option `--all-variants` will show these variants while the option `--best-variants` will omit them.

## Specific regions

By default, all genomic regions are shown with the option `--all-regions`. However, variants outside a gene (as annotated in the GFF file) can be removed from the output with the option `--gene-regions`. Variants are shown if their position in the VCF file are inside a gene (hence INDELs staring just before a gene will not be selected).

## <a name="scores"></a>Variant scores

The variant score displayed at the Score column of the CSV file is made of four different percentages separated by a ':' corresponding respectively to:

- the percentage of samples from group 1 (first name/identifier in CSV file or parents) for which the allele is a true variant
- the percentage of samples from group 1 (first name/identifier in CSV file or parents) for which the allele is a true reference
- the percentage of samples from group 2 (second name/identifier in CSV file or offspring) for which the allele is a true variant
- the percentage of samples from group 2 (second name/identifier in CSV file or offspring) for which the allele is a true reference

Note that when the `--comparison` option is not used, the third and fourth percentages are equal to the first and second.

As a result, fixed variant can either have the first and the third percentage equal to 0 (fixed reference in the population) or the second and the fourth equal to 0 (fixed variant in the population).

Differentially prevalent variants can either have the first and the fourth percentages different from 0, or the second and the third  percentages different from 0 (variant present in at least one sample from one group and absent from at least one sample from the other group).

Percentages within each group do not add up to 100 when there are samples that could not be genotyped (because of insufficient depth or a mixed variant).

# Additional information

## Up and downstream genomic sequences

The genomic sequences that are upstream or downstream of the reference allele can be shown in the Ref column of the output CSV file with respectively `--nucl-window-before` and `--nucl-window-after`.

## Overlapping features

Variant inside features are annotated given the feature type available in the third column of the GFF file. If the feature type contains the keyword 'gene', their annotation is saved in the Annotation column of the CSV file. If the feature type is 'CDS', the automatic translation will be processed and the CDS and amino acid sequences associated with the CDS annotation will be added in the CDSref, CDSalt, AAref and AAalt columns of the CSV file.
All the feature types detected are stored after the variant type in the Type column of the CSV file.

## Automatic translation

When the variant is inside a CDS, `varif` will predict the new CDS sequence affected by the variant by replacing the reference sequence by the mutation. However, knowing where to start the CDS can be difficult in the case of a mutation affecting the very first codon of the CDS. In that matter, several further modifications are added to produce the most biologically relevant CDS:

- When the mutated sequence is shorter than the reference sequence:
  - If the very first bases of the initial CDS are not recovered after adding the mutation, bases are added from the upstream sequence to lead to a mutated CDS of the same size than the initial one
  - Otherwise, bases from the downstream sequence are added to match the initial CDS size

- When the mutated sequence is longer than the reference sequence:
  - Bases located in the upstream region of the CDS are not included in the mutated CDS
  - All other bases from the mutation are kept (including those from the downstream region) and bases from the downstream sequence after the mutation are added to make the length of the mutated CDS a multiple of 3

Amino acids affected by the reference and alternate sequence are obtained by translating the initial and mutated CDS respectively. Additional codons can be added in a chosen window with `--prot-window-before` and `--prot-window-after` options; limited to the first and last codon (or when a stop codon is obtained). The position of the codon is then calculated from their rank from the beginning of the initial and mutated CDS and shown in the AAref and AAalt columns of the CSV file. The codon positions are associated with the total number of codons in each CDS, including the stop codon.

# Limitations

## Multi-variant CDS

As `varif` considers only one variant at a time, the wrong mutated CDS can be translated if two different mutations occur in the same one.

## CDS prediction accuracy

Since `varif` does not consider bases from a mutation if they reach the upstream sequence (before the first codon), a longer mutated CDS starting before the initial first codon could be missed.

Similarly, the end of the mutated CDS is not precisely determined but estimated by the addition of enough bases to reach at least the initial CDS size or a multiple of 3.

# Examples

1. Save all possible variants regardless of their VAFs.
    ```bash
    varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
        --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
        --fixed --all-variants --all-regions \
        --out-filename filtered-variants \
        --output-vcf
    ```

2. Save variants falling in a gene regardless of their VAFs with 10 bases upstream and downstream of the variants.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --fixed --all-variants --gene-regions \
       --nucl-window-before 10 --nucl-window-after 10 \
       --out-filename filtered-variants
   ```
   
3. Save variants falling in a gene and differentially present (at least one true reference and true variant in all samples) with protein sequences including 5 amino acids before and after the variant.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --no-fixed --best-variants --gene-regions \
       --prot-window-before 5 --prot-window-after 5 \
       --out-filename filtered-variants
   ```
   
4. Save variants falling in a gene and differentially present between families (at least one true reference in a family and a true variant in the other) with protein sequences including 5 amino acids before and after each variant.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --no-fixed --best-variants --gene-regions \
       --ped my_ped.ped --comparison families \
       --prot-window-before 5 --prot-window-after 5 \
       --out-filename filtered-variants
   ```
   
4. Save variants falling in a gene and differentially present in lineages (at least one true reference in a parent and a true variant in an offspring and inversely) only if there are at least 80% of non missing genotypes in both groups.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --no-fixed --best-variants --gene-regions \
       --ped my_ped.ped --comparison lineages \
       --max-missing 0.2 \
       --prot-window-before 5 --prot-window-after 5 \
       --out-filename filtered-variants
   ```
