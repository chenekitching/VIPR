# VIPR (Variant Interpretation and Prioritisation using Random forest) 

## Nextflow pipeline to filter, annotate, and prioritise genetic variants

## Description
The pipeline takes a gzipped VCF file as input, splits it by chromosome, and performs the filtering annotation, and prioritisation steps steps in parallel. The prioritisation step uses a random forest machine learning model to prioritise the variants according to the ACMG criteria.

## Dependancies
1. [Nextflow](https://www.nextflow.io)
2. [Annovar](http://annovar.openbioinformatics.org/en/latest/)
   - The following databases need to be installed:
        - refGene
        - ensGene
        - avsnp150
        - gnomad211_genome (for hg19) or gnomad40_genome (for hg38)
        - dbnsfp42a
        - clinvar_20221231
        - intervar_20180118 
3. [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)

## Usage
Parameters can either be specified in the nextflow.config file, or directly in the command line as per the example below.
Simple use case example:
```bash
nextflow run vipr.nf --vcf_file path/to/vcf --indx_file path/to/index --build hg19 --annovar_db path/to/humandb --outdir path/to/output/directory
```
## Parameters
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| `--vcf_file`    |    | bgzipped VCF file |
| `--indx_file`     |   | csi index for vcf.gz file |
| `--annovar_db`    |    |  Folder with annovar databases |
| `--build`    |  hg19 |  Reference genome |
| `--outdir` | results  |  Output folder |
| `--mod` | $projectDir/weighted_rf_18-07.rds  |  Location of the machine learning model used for variant prioritisation |

## Output
  | Type      | Description     |
  |-----------|---------------|
  | merge.txt      | Unsorted variants with ML prediction appended to annovar output |
  | prioritised_file.txt  | Sorted variants with ML prediction appended to annovar output |
  | trimmed_prioritised_file.txt  | Variants from prioritised_file.txt with > 0.5 probability of pathogenicity |
