# vcf_annotate

## Description
A prototype VCF parser and ExAC annotation tool. In the command line interface, the user supplies a VCF file (v4.1) to the module, which returns an annotated CSV file, with additional variant information columns extracted via the ExAC database API.

ExAC API can be found here:
>  http://exac.hms.harvard.edu/

Each variant has the following information: 
- *var_type*: type of variation (CNV, indel, etc.)
- *effect*:  amino acid effect (missense, silent, etc.) 
- *depth*: depth of sequence coverage at the site of variation
- *num_reads*: number of reads supporting the variant versus those supporting reference reads
- *exac_allele_freq*: allele frequency, from ExAC

Note: If there are multiple variant effects, the annotation is the *most deleterious* one. Script contains the an ordered list of variant effect severity, sourced from Ensembl.

## Usage
Command line:
```
python3 vcf_annotate.py —INPUT.vcf —OUTPUT_FILE_NAME
```


## Dependencies
See requirements.txt
