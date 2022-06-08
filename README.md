# varif
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
  **Used to annotate the variants**
- *CDS*
  **Used to translate sequences affected by the variants**

## FASTA (`-fasta`)

The FASTA file containing one sequence per chromosome.

# Process

Use `varif` to filter and rank variants submitted through a VCF file.

For each variant in each sample, `varif` will calculate the Variant Alternate Frequency (VAF) and classify it as a true variant or not based on filters applied with the available options. The values of VAFs for each variant are compared to find sites likely to be differentially distributed among samples.

# Output

Variants passing the filters will be sent to the standard output with the option `--show` and/or to a VCF file if the option `--filtered-vcf` is set. Additionally, a CSV file containing each individual variant (variants from the same chromosome location are separated) is generated if the option `--filtered-csv` is set.

## CSV content

The output CSV file separated by semicolons contains:

- Chromosome: The chromosome name
- Position: The starting position of the variant
- Type: The type of variant (SNP, INDEL...)
- Ref: The reference sequence
- Alt: The alternate sequence
- CDSref: The windowed CDS reference sequence including the variant followed by the CDS identifier\*
- CDSalt: The windowed CDS alternate sequences including the variant\*
- AAref: The amino acid sequence resulting from the windowed reference sequence\*
- AAalt: The amino acid sequence resulting from the windowed alternate sequences\*
- Annotation: The annotation described in the GFF file followed by the gene identifier\*\*
- Score: The score of the variant (see [Process](#process))
- Sample#1: Proportion of alternate allele for sample 1
- ...
- Sample#n: Proportion of alternate allele for sample n
- Con: Sample(s) having the reference sequence
- Mut: Sample(s) having the alternate sequence
- Mix: Sample(s) having a mixture of sequences
- Und: Sample(s) having not enough depth at this position

\*: Applies for CDS regions, otherwise *NA*
\*\*: Applies for CDS and intronic regions, otherwise *NA*

# Filters

## Minimal depth

The variant allele frequencies are calculated only if the sample allele depth (the sum of all the reads REF and ALTs) is equal or superior to the value provided with the option `--depth`, with a default value of 5.

## Cut-off proportions

There are 2 cut-offs of minimal Alternate Allele Frequency (minAAF with `ratio-alt`) and maximal Reference Allele Frequency (maxRAF with `ratio-no-alt`) used by `varif` to determine if variants are differentially expressed in the population. By default, minAAF is equal to 0.8 and minRAF to 0.2.

For each variant at a chromosomal location, the variant allele frequency (VAF) is calculated using the different allele depths (AD):

***VAF*** = (*variant AD*) / (*other variants AD* + *reference AD*)

If the VAF is above minAAF, it should then be considered a true variant, while if the VAF is below minRAF, it should be considered a true reference.

The combination of VAFs are used to classify variants that are:

- Differentially distributed among samples:
  At least one VAF of the samples is below or equal to maxRAF and at least one VAF is above or equal to minAAF.
- Fixed in the population:
  All VAFs of the samples are above or equal to minAAF.
- Not differentially distributed among samples, for all other cases.

## Specific regions

By default, all genomic regions are shown with the option `--all-regions`. However, variants outside a gene (as annotated in the GFF file) can be removed from the output with the option `--gene-regions`. Variants are shown if their position in the VCF file are inside a gene (hence INDELs staring just before a gene will not be selected).

## Variant scores

Variants can also be filtered by their score which are calculated differently whether:

- The variant is differentially expressed among samples:
  *Score* = *maxVAF* / *minVAF*
  If the minVAF is equal to 0, a maximum can be set with the option `--max-score` (99999 by default)
- The variant is fixed in the population:
  *Score* = 0​
  All these variants can be omitted from the output by adding the option `--no-fixed`
- All other cases:
  *Score* = -*maxVAF* / *minVAF*
  Similarly, if the minVAF is equal to 0, a minimum can be set with the option `--min-score` (-99999 by default)

# Additional information

## Up and downstream genomic sequences

The genomic sequences that are upstream or downstream of the reference allele can be shown in the Ref column of the output CSV file with respectively `--nucl-window-before` and `--nucl-window-after` .

## Multiple features

Variant inside genes are annotated given the gene annotation available in the GFF file. Several gene annotations are separated by a ':' , as well as the associated CDS and amino acid sequences when the variant is inside a CDS.

## Automatic translation

When the variant is inside a CDS (i. e. its first position is inside a CDS), this feature will predict the protein sequence affected by the whole variant plus bases before and after included in a chosen window with `--prot-window-before` and `--prot-window-after` options.
Windows are expressed in amino acid before the and after the whole variant in the direction of the translation.
The translation will stop if:

- The window is completely translated

- The end of the CDS is reached

- A stop codon is obtained

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

3. Save variants falling in a gene and differentially expressed (positive VAF) with protein sequences including 5 amino acids before and after the variant.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       --depth 5 --ratio-alt 0.8 --ratio-no-alt 0.2 \
       --no-fixed --best-variants --gene-regions \
       --prot-window-before 5 --prot-window-after 5 \
       --filtered-csv filtered-variants.csv \
       --filtered-vcf filtered-variants.vcf
   ```

   

