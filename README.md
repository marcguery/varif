# varif
A filtering and annotating program for VCF-formatted variants.

# Download and installation

```bash
git clone https://github.com/marcguery/varif.git
cd varif
pip3 install .
varif -h
```

# Inputs

All input files must have the same chromosome names.

## VCF

A VCF file formatted as described in `samtools` specifications at [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

Mandatory metadata headers are:

-  #CHROM
- POS
- ID
- REF
- ALT
- QUAL
- FILTER
- INFO

## GFF

A GFF3 file formatted as described in `ensembl` specifications at [https://www.ensembl.org/info/website/upload/gff3.html](https://www.ensembl.org/info/website/upload/gff3.html).

## FASTA

The FASTA file containing one sequence per chromosome.

# Process

Use `varif` to filter and rank variants submitted through a VCF file.

For each variant in each sample, `varif` will calculate the proportion of alternate reads. 
Then, it will calculate a score based on the ratio of the proportion of alternate reads of the most variant sample above that of the less variant sample. This score will be then modified considering chosen values of ratio to call a variant true positive or true negative (see [Cut-off proportions](#cut\\-off-proportions)).
As a result, the highest scores will be associated with sites likely to be differentially distributed among samples.

## Output

These scores will be saved in a CSV file whose headers are:

- Chromosome: The chromosome name
- Position: The starting position of the variant
- Type: The type of variant (SNP, INDEL...)
- Ref: The reference sequence
- Alt: The alternate sequence
- AAref: The amino acid resulting from the reference sequence followed by the CDS identifier
- AAalt: The amino acid resulting from the alternate sequence
- Annotation: The annotation described in the GFF file followed by the gene identifier
- Score: The score of the variant (see [Process](#process))
- Sample#1: Proportion of alternate allele for sample 1
- ...
- Sample#n: Proportion of alternate allele for sample n
- Con: Samples having the reference sequence
- Mut: Sample having the alternate sequence
- Mix: Samples having a mixture of sequences
- Und: Samples having not enough depth at this position

# Filters

## Minimal depth

The depth below which all the read counts of the corresponding sample are set to null.
For all these samples, proportions are NA values.

## Cut-off proportions

There are 2 cut-offs used by `varif`, the first one (**cut-off #1**) being the minimal proportion to be considered a true variant, the second one (**cut-off #2**) the maximal proportion to be considered a true reference.

These cut-offs are used to modify the score in order to know whether the site is likely to be:

- Differentially distributed among samples, hence the less variant is below **cut-off #2** and the most variant is above **cut-off #1**. *The score is not modified*.
- Fixed in the population, hence the less variant site is still above **cut-off #1**. *The score is set to 0*.
- Not differentially distributed among samples, for all other cases. *The opposite of the score becomes the score*.

## Specific regions

Variants can be filtered by their location, i. e. either inside or outside a gene as annotated in the GFF file. Variants are shown if their position in the VCF file are inside a gene (hence INDELs staring just before a gene will not be selected).

## Variant types

Variants can also be filtered by their score; you can show all variants whose score is:

- Any value: all variants
- Equal or superior to 0: fixed and differentially expressed variants
- Superior to 0: differentially expressed variants.

# Additional information

## Annotation

Variant inside genes are annotated given the gene annotation available in the GFF file. 

Variants outside genes will have a *NA* instead.

## Automatic translation

When the variant is inside a CDS (i. e. its first position is inside a CDS), this feature will predict the protein sequence affected by the variant.
Note that this prediction could not be reliable when:

- The start of the variant is located 2 bases or less before the end of the CDS
- The end of the variant is located 2 bases or less before the end of the CDS or is longer than it