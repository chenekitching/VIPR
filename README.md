# nf-variant-prioritisation

## Nextflow pipeline to filter, annotate, and prioritise genetic variants

## Description
The pipeline takes a gzipped VCF file as input, splits it by chromosome, and performs the filtering and annotation steps in parallel. The final process uses a random forest machine learning model to prioritise the variants according to the ACMG criteria.

## Dependancies
1. Nextflow
2. Annovar
3. Singularity

## Usage
