# nf-variant-prioritisation

## Nextflow pipeline to filter, annotate, and prioritise genetic variants

## Description
The pipeline takes a gzipped VCF file as input, splits it by chromosome, and performs the filtering and annotation steps in parallel. The final process uses a random forest machine learning model to prioritise the variants according to the ACMG criteria.

## Dependancies
1. Nextflow
2. Annovar
3. Singularity

## Usage
Simple use case example:
```bash
nextflow run vp.nf --vcf_file path/to/vcf --indx_file path/to/index --build hg19 --annovar_db path/to/humandb --outdir path/to/output/directory
```
## Parameters
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| `--vcf_file`    |  ``  | bgzipped VCF file |
| `--indx_file`     |  `` | csi index for vcf.gz file |
| `--annovar_db`    |  ``  |  Folder with annovar databases |
| `--build`    |  hg19 |  Reference genome |
| `--outdir` | results  |  Output folder |

## Output
  | Type      | Description     |
  |-----------|---------------|
  | merge.txt      | Unsorted variants with ML prediction appended to annovar output |
  | prioritised_file.txt  | Sorted variants with ML prediction appended to annovar output |
  | trimmed_prioritised_file.txt  | Variants from prioritised_file.txt with > 0.5 probability of pathogenicity |
