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

# Overview

Use `varif` to filter and annotate variants submitted through a VCF file.

For each allele in each sample, `varif` will calculate the Allele Sample Proportion (ASP) and classify it as a true mutation or not based on filters applied with the available options. The values of all ASPs are compared to find sites likely to be differentially mutated between group of samples.

Alleles passing the filters will be stored in a CSV file (multiple alleles at the same position are separated) generated with the name (without extension) set by the option `-outfilename`. If the option `--output-vcf` is set, the corresponding filtered variants (multiple alleles at the same position are merged) from the VCF file will be saved as well. Each combination of  groups compared will lead to a unique CSV (and VCF) file (see [Sample groups](#groups)).

Here is an example of the  `varif` main output obtained with default parameters:

| Chromosome | Position | Type                     | Ref  | Alt  | CDSref           | CDSalt   | AAref      | AAalt       | Annotation                   | Proportions       | S1       | S2       | S3       |
| ---------- | -------- | ------------------------ | ---- | ---- | ---------------- | -------- | ---------- | ----------- | ---------------------------- | ----------------- | -------- | -------- | -------- |
| Chr1       | 123123   | SNP                      | C    | G    | NA               | NA       | NA         | NA          | NA                           | `067:033:067:033` | 0.912643 | 0.000000 | 1.000000 |
| Chr1       | 456456   | SNP,gene,mRNA            | T    | A    | NA               | NA       | NA         | NA          | Unknown function (gene_id_1) | `033:033:033:033` | 0.540000 | 0.000000 | 1.000000 |
| Chr2       | 123123   | SNP,CDS,exon,gene,mRNA   | G    | T    | AAGCG (CDS_id_2) | AATCG    | L (45/458) | I (45/458)  | Unknown function (gene_id_2) | `000:067:000:067` | NA       | 0.000000 | 0.140000 |
| Chr3       | 123456   | INDEL,CDS,exon,gene,mRNA | G    | GATA | ATGAT (CDS_id_3) | ATGATAAT | D (52/123) | DN (52/124) | Unknown function (gene_id_3) | `000:100:000:100` | 0.000000 | 0.120000 | 0.000000 |

## CSV content

The output CSV file separated by semicolons contains:

- ***Chromosome*** : The chromosome name
- ***Position*** : The starting position of the allele
- ***Type*** : The type of allele (SNP, INDEL...) followed by genomic features found
- ***Ref*** : The reference sequence +/- up and downstream sequences (separated by '|')
- ***Alt*** : The alternate sequence
- ***CDSref*** : The windowed CDS reference sequence including the allele followed by the CDS identifier\*
- ***CDSalt*** : The windowed CDS alternate sequences including the allele\*
- ***AAref*** : The amino acids resulting from the windowed reference sequence followed by the codon position\*
- ***AAalt*** : The amino acids resulting from the windowed alternate sequences followed by the codon position\*
- ***Annotation*** : The annotation described in the GFF file followed by the gene identifier\*\*
- ***Proportions*** : The proportion of mutated and reference alleles in the population (see [Allele Population Proportions](#app))
- The remaining columns contain the allele proportion of each sample (ASPs)

\*: Applies for CDS regions, otherwise *NA*

\*\*: Applies for transcribed regions, otherwise *NA*

Note that when the CDSalt field is empty and that the AAalt field contains '*(0/0)*', this means that both reference and mutated CDS sequences are identical (this can happen when an INDEL occurs at the edge of a CDS and an intron).

## <a name="groups"></a>Sample groups

Samples groups are used to determine which allele is differentially distributed based on the presence of true mutated samples and true reference samples between two groups. It is possible to separate samples into groups and look for alleles differentially distributed between those groups by providing the option `--comparison` with one of the keywords *families*,*selfself*, *lineages*, or *all* (the default). Samples are grouped into families or parental-offspring relationships according to the metadata from the PED file. Allele differential distribution can be compared:

- between every combination of families (*families*),
- within all members of each family (*selfself*),
- between direct parents and offspring  (*lineages*),
- or all-against-all VCF samples (*all*).

The relationship between the two groups is symmetrical, meaning that alleles present in either group are candidates for being considered differentially distributed (see [Allelic distribution](#alleledistrib) for more information).

# Filters

## Allele Sample Proportions (ASP)

There are 2 cut-offs of minimal alternate ASP (minVSP with `--ratio-alt`) and maximal reference ASP (maxRSP with `--ratio-ref`) used by `varif` to determine if alleles are differentially distributed in a population. By default, minVSP is equal to 0.8 and maxRSP to 0.2.

For each allele and for each sample, the ASP is calculated using the different allele read depths (ARD):

***ASP*** = (*alternate ARD*) / (*all alternates ARD* + *reference ARD*)

If the ASP is above minVSP, the locus is considered to carry the alternate allele only, while if the ASP is below maxRSP, it is considered to not carry the alternate allele at all. When the ASP is between maxRSP and minVSP, the genotype cannot be called and the locus is mixed.

## Minimal depth

The ASPs are calculated only if the sample read depth (the sum of all the reads REF and ALTs) is equal or superior to the value provided with the option `--depth`, with a default value of 5.

## Missing data

Variants with too many missing genotypes (because of insufficient depth or a mixed ASP) can be filtered out with `--max-missing` (default 1) which discards alleles if one group has more missing genotypes than the selected proportion.

## Rare alleles

Alleles that are too rare can be filtered out by the option `--min-mutated` (default 0) which will discard alleles that are present in less than this proportion of the number of samples in both group.

## Group similarity

Group-specific alleles can be selected with the option `--max-similarity` (default 1) by providing the maximal ratio of the mutant proportion of the less mutated group over the mutant proportion of the most mutated group. 

For example, a value of 0.2 will select only alleles with a proportion of mutated samples that is at least 5 times higher in one group compared to the other.

## Genetic regions

By default, all genomic regions are included with the option `--all-regions`. However, alleles outside a gene (all of the reference bases are located outside of a gene or in an intron) can be removed from the output with the option `--gene-regions`.

# <a name="app"></a>Allele Population Proportions

The combination of ASPs are used to classify alleles that are:

- Differentially distributed between two groups:
  At least one ASP of one group is below or equal to maxRSP (does not contain the alternate allele) and at least one ASP of the other group is above or equal to minVSP (contains the alternate allele).
- Fixed in the population:
  All ASPs of the samples are above or equal to minVSP or below or equal to maxRSP. Either the option `--all-variants` or `--fixed` will show these alleles.
- Not differentially distributed among samples, for all other cases. The option `--all-variants` will show these alleles while the option `--best-variants` will omit them. Note that the whole variant (with all alternate alleles merged in a single line) will be saved in the VCF file even if only one allele is passing the filters.

The exact proportion of mutated and reference alleles in both groups, refer to as 'Allele Population Proportions', is displayed at the Proportions column of the CSV file. There are 4 different percentages separated by a ':' corresponding respectively to:

- the percentage of samples from group 1 (first name/identifier in CSV file or parents) for which the locus *is* the alternate allele
- the percentage of samples from group 1 (first name/identifier in CSV file or parents) for which the locus *is not* the alternate allele
- the percentage of samples from group 2 (second name/identifier in CSV file or offspring) for which the locus *is* the alternate allele
- the percentage of samples from group 2 (second name/identifier in CSV file or offspring) for which the locus *is not* the alternate allele

Note that when the `--comparison` option is set to *all*, the third and fourth percentages are equal to the first and second.

As a result, a fixed allele can either have the first and the third percentage equal to 0 (fixed reference in the population) or the second and the fourth equal to 0 (fixed allele in the population).

Differentially distributed alleles can either have the first and the fourth percentages different from 0, or the second and the third  percentages different from 0 (allele present in at least one sample from one group and absent from at least one sample from the other group).

Percentages within each group do not add up to 100 when there are samples that could not be genotyped (because of insufficient depth or a mixed allele).

# Sequence annotations

## Up and downstream genomic sequences

The genomic sequences that are upstream or downstream of the reference allele can be shown in the Ref column of the output CSV file with respectively `--nucl-window-before` and `--nucl-window-after`.

## Overlapping features

Alleles inside features are annotated given the feature type available in the third column of the GFF file. If the feature type contains the keyword 'gene', their annotation is saved in the Annotation column of the CSV file. If the feature type is 'CDS', the automatic translation will be processed and the CDS and amino acid sequences associated with the CDS annotation will be added in the CDSref, CDSalt, AAref and AAalt columns of the CSV file.
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

Amino acids affected by the reference and alternate sequence are obtained by translating the initial and mutated CDS respectively. Additional codons can be added in a chosen window with `--prot-window-before` and `--prot-window-after` options; limited to the first and last codon (or when a stop codon is obtained). The position of the codon is then calculated from their rank from the beginning of the initial and mutated CDS and shown in the AAref and AAalt columns of the CSV file. The codon positions are associated with the total number of codons in each CDS, including the stop codon.

# Multiprocessing and variant batches

The most time consuming step of the pipeline is when the ASPs are calculated. To save on memory, variants are processed by chunk whose number is determined from an upper limit of variants processed by chunk set by `--chunk-size` (5000 variants by chunk by default). Several chunks can be processed at the same time by increasing the number of cores with `--ncores` (1 by default). The number of chunks is increased to match the number of cores if necessary.

# Limitations

## Multi allelic CDS

As `varif` considers only one allele at a time, the wrong mutated CDS can be translated if two different mutations occur in the same one.

## CDS prediction accuracy

Since `varif` does not consider bases from a mutation if they reach the upstream sequence (before the first codon), a longer mutated CDS starting before the initial first codon could be missed.

Similarly, the end of the mutated CDS is not precisely determined but estimated by the addition of enough bases to reach at least the initial CDS size or a multiple of 3.

## Mutation overlapping CDS and introns

if an INDEL includes the edge between a CDS and an intron, the resulting protein could include false codons from the intronic regions as for now no mutated base within an intronic region are removed.

# Examples

1. Save all possible alleles regardless of their ASPs, processing in parallel 4 batches of 1000 variants.
    ```bash
    varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
        -outfilename ./out/filtered-variants \
        --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
        --fixed --all-variants --all-regions \
        --output-vcf --chunk-size 1000 --ncores 4
    ```

2. Save alleles falling in a gene regardless of their ASPs with 10 bases upstream and downstream of the alleles.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --fixed --all-variants --gene-regions \
       --nucl-window-before 10 --nucl-window-after 10
   ```
   
3. Save alleles falling in a gene and differentially distributed (at least one sample with the allele and one sample without it among all samples) with protein sequences including 5 amino acids before and after the allele.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --no-fixed --best-variants --gene-regions \
       --prot-window-before 5 --prot-window-after 5
   ```
   
4. Save alleles falling in a gene and differentially distributed between families (at least one sample with the allele in a family and one sample without it in another one) with protein sequences including 5 amino acids before and after each allele.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --no-fixed --best-variants --gene-regions \
       --ped my_ped.ped --comparison families \
       --prot-window-before 5 --prot-window-after 5
   ```
   
4. Save alleles falling in a gene and differentially distributed in lineages (at least one sample with the allele in a parent and one sample without it in an offspring and inversely) only if there are at least 80% of non missing genotypes in both groups.

   ```bash
   varif -vcf my_vcf.vcf -gff my_gff.gff -fasta my_fasta.fasta \
       -outfilename ./out/filtered-variants \
       --depth 5 --ratio-alt 0.8 --ratio-ref 0.2 \
       --no-fixed --best-variants --gene-regions \
       --ped my_ped.ped --comparison lineages \
       --max-missing 0.2 \
       --prot-window-before 5 --prot-window-after 5
   ```
