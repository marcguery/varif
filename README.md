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

Mandatory metadata headers are :

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
Then, it will calculate a score based on the ratio of the proportion of alternate reads of the most variant sample above that of the less variant sample.
As a result, the highest scores will be associated with sites probably differentially distributed among samples.

These scores will be saved in a CSV file. Subset the results by either the score values (showing fixed variants, differentially expressed or all variants) or the region where they were found (inside or outside genes).

# Filters

## Minimal depth

The depth below which all the read counts of the corresponding sample are set to null.
For all these samples, proportions are NA values.

## Cut-off proportions

There are 2 cut-offs used by `varif`, the first one (**cut-off #1**) being the minimal proportion to be considered a true variant, the second one (**cut-off #2**) the maximal proportion to be considered a true reference.

These cut-offs are used to modify the score in order to know whether the site is likely to be:

- Differentially distributed among samples, hence the less variant is below **cut-off #2** and the most variant is above **cut-off #1**. *The score is not modify*.
- Fixed in the population, hence the less variant site is still above **cut-off #1**. *The score is set to 0*.
- Not differentially distributed among samples, for all other cases. *The opposite of the score becomes the score*.



