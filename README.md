# VIPR (Variant Interpretation and Prioritisation using Random forest) 

## Nextflow pipeline to filter, annotate, and prioritise genetic variants

## Description
The pipeline takes a gzipped VCF file as input, splits it by chromosome, and performs the filtering annotation, and prioritisation steps steps in parallel. The prioritisation step uses a random forest machine learning model to prioritise the variants according to the ACMG criteria.

## Dependencies
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
### Example 
To run the example, simply execute the following command, and specify the path to your ANNOVAR database (either in the nextflow.config file or as --annovar_db path/to/humandb in the command below):
```bash
nextflow run vipr.nf
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

## Shiny UI
To interactively explore the results in the prioritised_file.txt, the Shine user interface can be launched by running the following command. The file name and genome build are specified as parameters. The Shiny app is executed within a Singularity container, so there is no need to install any packages.

```bash
singularity exec --bind /host/path:/container/path shiny_cont.sif Rscript -e 'shiny::runApp("/mnt/vp_4.R", launch.browser = TRUE)'
```
Here, /host/path is the absolute path on your host system that contains the shiny script (vipr_shiny.R). /container/path is the directory inside the container where the host path will be mounted. To confirm the /container/path, you can execute the code below:
To enter the container, run
```bash
singularity shell shiny_cont.sif
```
To see what the /container/path is, sun
```bash
pwd
```
