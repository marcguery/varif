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

# Process

Use `varif` to filter and rank variants submitted through a VCF file.

For each variant in each sample, `varif` will calculate the Variant Alternate Frequency (VAF) and classify it as a true variant or not based on filters applied with the available options. The values of VAFs for each variant are compared to find sites likely to be differentially distributed among samples.

# Output

Variants passing the filters will be sent to the standard output by default or to a VCF file if the option `--filtered-vcf` is set. Additionally, a CSV file containing each individual variant (variants from the same chromosome location are separated) is generated if the option `--filtered-csv` is set.

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

## Control groups

It is possible to separate samples into 2 groups (positive and negative control groups) that are used to determine which variant is differentially expressed between the 2 groups based on the allele frequencies. By default, all samples belong to the same group which implies that the positive and negative control groups are the same. The option `--control` can be used to specify the negative control group with a comma separated list of sample identifiers corresponding to their order in the VCF header, starting by 1. The relationship between the 2 groups is symmetrical, meaning that variants detected in both groups are candidates for being considered differentially expressed.

## Allele frequencies

There are 2 cut-offs of minimal Alternate Allele Frequency (minAAF with `--ratio-alt`) and maximal Reference Allele Frequency (maxRAF with `--ratio-no-alt`) used by `varif` to determine if variants are differentially expressed in the population. By default, minAAF is equal to 0.8 and maxRAF to 0.2.

For each variant at a chromosomal location and for each sample, the variant allele frequency (VAF) is calculated using the different allele depths (AD):

***VAF*** = (*variant AD*) / (*other variants AD* + *reference AD*)

If the VAF is above minAAF, it should then be considered a true variant, while if the VAF is below maxRAF, it should be considered a true reference.

The combination of VAFs are used to classify variants that are:

- Differentially distributed between positive and negative control samples:
  At least one VAF of the **negative** control samples is a true reference (below or equal to maxRAF) and at least one VAF of the **positive** control samples is a true variant (above or equal to minAAF).
  The rule applies symmetrically, when at least one VAF of the **positive** control samples is a true reference and at least one VAF of the **negative** control samples is a true variant.
- Fixed in the population:
  All VAFs of the samples are above or equal to minAAF or below or equal to maxRAF. Either the option `--all-variants` or `--fixed` will show these variants.
- Not differentially distributed among samples, for all other cases. The option `--all-variants` will show these variants while the option `--best-variants` will omit them.

## Specific regions

By default, all genomic regions are shown with the option `--all-regions`. However, variants outside a gene (as annotated in the GFF file) can be removed from the output with the option `--gene-regions`. Variants are shown if their position in the VCF file are inside a gene (hence INDELs staring just before a gene will not be selected).

## <a name="scores"></a>Variant scores

The variant score displayed at the Score column of the CSV file is made of four different percentages separated by a ':' corresponding respectively to:

- the percentage of positive control samples for which the allele is a true variant
- the percentage of positive control samples for which the allele is a true reference
- the percentage of negative control samples for which the allele is a true variant
- the percentage of negative control samples for which the allele is a true reference

As a result fixed variant can be either have the first and the third percentage equal to 0 (fixed reference in the population) or the second and the fourth equal to 0 (fixed variant in the population).

Differentially expressed variants can either have the first and the fourth percentages different from 0 (variant present in at least one sample from the positive control group and absent from at least one sample from the negative control group) or the second and the third (variant present in at least one sample from the negative control group and absent from at least one sample from the positive control group).

Percentages within each group do not add up to 100 when there are samples that have not sufficient read depth to calculate an allele frequency.

# Additional information

## Up and downstream genomic sequences

The genomic sequences that are upstream or downstream of the reference allele can be shown in the Ref column of the output CSV file with respectively `--nucl-window-before` and `--nucl-window-after` .

## Overlapping features

Variant inside features are annotated given the feature type available in the third column of the GFF file. If the feature type contains the keyword 'gene', their annotation is saved in the Annotation column of the CSV file. If the feature type is 'CDS', the automatic translation will be realised and the CDS and amino acid sequences associated with the CDS annotation will be added in the CDSref, CDSalt, AAref and AAalt columns of the CSV file.
All the feature types detected are stored after the variant type in the Type column of the CSV file.

## Automatic translation

When the variant is inside a CDS, `varif` will predict the new CDS sequence affected by the variant by replacing the reference sequence by the mutation. However, knowing where to start the CDS can be difficult in the case of a mutation affecting the very first codon of the CDS. In that matter, several further modifications are added to produce the most biologically relevant CDS. Bases from the upstream sequence (before the first codon) are added to complement the mutated CDS it it lacks the very first bases so that the initial and mutated have the same size and their resulting protein can be easilly compared. Similarly, bases from the mutation belonging to the upstream sequence are removed from the mutated CDS to match the initial CDS size.

More generally:

-  When the mutation is shorter than the reference sequence:
  - If the very first bases of the initial CDS are not recovered after adding the mutation, bases are added from the upstream sequence to lead to a mutated CDS of the same size than the initial one.
  - Otherwise, bases from the downstream sequence are added to match the initial CDS size
- When the mutation is longer than the reference sequence:
  - Bases located in the upstream region of the CDS are removed
  - All other bases from the mutation are kept (including those from the downstream region) and bases from the downstream sequence after the mutation are added to make the length of the mutated CDS a multiple of 3

Amino acids affeced by the reference and alternate sequence are obtained by translating the initial and mutated CDS respectively. Additional codons can be added in a chosen window with `--prot-window-before` and `--prot-window-after` options; limited to the first and last codon (or when a stop codon is obtained). The position of the codon is then calculated from their rank from the beginning of the initial and mutated CDS and shown in the AAref and AAalt columns of the CSV file. The codon positions are associated with the total number of codons in each CDS, including the stop codon.

# Limitations

## Multi-variant sample

As `varif` considers only one variant at a time for each sample which can lead to the wrong mutated CDS being translated.

## CDS prediction

Since `varif` does not consider bases from a mutation if they reach the upstream sequence (before the first codon), a longer mutated CDS starting before the initial first codon could be missed.

Similarly, the end of the mutated CDS is not precisely determined but estimated by the addition of enough bases to reach at least the initial CDS size or a multiple of 3.

# Examples

1. Save all possible variants regardless of their VAFs.
    ```bash
    varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
        --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
        --fixed --all-variants --all-regions \
        --filtered-csv filtered-variants.csv \
        --filtered-vcf filtered-variants.vcf
    ```

2. Save variants falling in a gene regardless of their VAFs with 10 bases upstream and downstream of the variants.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --fixed --all-variants --gene-regions \
       --nucl-window-before 10 --nucl-window-after 10 \
       --filtered-csv filtered-variants.csv \
       --filtered-vcf filtered-variants.vcf
   ```

3. Save variants falling in a gene and differentially expressed (all samples vs all other samples) with protein sequences including 5 amino acids before and after the variant.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --no-fixed --best-variants --gene-regions \
       --prot-window-before 5 --prot-window-after 5 \
       --filtered-csv filtered-variants.csv \
       --filtered-vcf filtered-variants.vcf
   ```

4. Save variants falling in a gene and differentially expressed (3rd & 4th samples in VCF file header vs all other samples) with protein sequences including 5 amino acids before and after the variant.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --no-fixed --best-variants --gene-regions \
       --control 3,4 \
       --prot-window-before 5 --prot-window-after 5 \
       --filtered-csv filtered-variants.csv \
       --filtered-vcf filtered-variants.vcf
   ```

