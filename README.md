# vcf_annotate

## Description
A prototype VCF parser and ExAC annotation tool. In the command line interface, the user supplies a VCF file (v4.1) to the module, which returns an annotated CSV file, with additional variant information columns extracted via the ExAC database API.

ExAC API can be found here:
>  http://exac.hms.harvard.edu/

Each variant has the following information: 
- *var_type*: type of variation (complex, insertion, deletion, etc.)
- *effect*:  genotype consequence (missense, synonymous, etc.) 
- *depth*: depth of sequence coverage at the site of variation
- *num_reads*: number of reads supporting the variant versus those supporting reference reads
- *percent_variant_support*: percent of reads supporting variant
- *exac_allele_freq*: allele frequency, from ExAC

Note: If there are multiple variant effects, the annotation is the *most deleterious* one. Script contains an ordered list of variant effect severity, sourced from Ensembl.

## Usage
Command line:
- Names in < > are for the user to replace with the appropriate file names.
```
python3 vcf_parser_prototype.py <INPUT.vcf> <OUTPUT_FILE_NAME>
```


## Dependencies
See environment.yml
