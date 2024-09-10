#!/usr/bin/env nextflow

/*conditional gnomad params: > gnomad 2 NA to hg19
    default version for hg38 is gnomad 4
*/
if ("${params.build}" == "hg19" ) {
    println "Using hg19."
    params.gnomad = "gnomad211_genome"
    params.af = "AF"
}
else if ("${params.build}" == "hg38" ) {
    params.gnomad = "gnomad40_genome"
    params.af = "gnomad40_genome_AF"
}


 log.info """\
    V I P R   P I P E L I N E
    =====================
    vcf          : ${params.vcf_file}
    index        : ${params.indx_file}
    outdir       : ${params.outdir}
    build        : ${params.build}
    """
    .stripIndent(true)

process SPLIT{

    input:
    path vcf
    path indx

    output:
    path 'chr*' 

    script:
    """
    bash $projectDir/split.sh $vcf $indx chr
    
    """
}

process FILTER {

    input:
    path split_vcf

    output:
    path "${split_vcf.baseName}.filt.vcf"


    """
    bcftools view -f PASS $split_vcf > ${split_vcf.baseName}.filt.vcf
    
    """

}

process ANNOTATE {

  input:
  path filtered_file 

  output:
  path "${filtered_file.baseName}.*"
  
  script:
  """
  table_annovar.pl $filtered_file ${params.annovar_db} -buildver ${params.build} -out ${filtered_file.baseName}.* -remove -protocol refGene,ensGene,avsnp150,${params.gnomad},dbnsfp42a,clinvar_20221231,intervar_20180118 -operation g,g,f,f,f,f,f -arg '-hgvs',,,,,, -nastring . -vcfinput -polish
  
  """
}

process ML {

  input:
  path mod
  path anno_file
  
  output:
  path 'prioritised_file.txt'
  
  script:
  """
  #!/usr/bin/env Rscript

# Load libraries ----------------------------------------------------------

library(dplyr)
library(readr)
library(modelr)
library(tidymodels)

# CommandArgs -------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

# Read files --------------------------------------------------------------
mod <- readRDS("${mod}")
input <- read_tsv("${anno_file}") %>%
  rowid_to_column("ID") %>%
  mutate(varID = paste(
    Chr,
    Start,
    End,
    Ref,
    Alt,
    sep = "_"
  ))

# Functions ---------------------------------------------------------------

get_annotations <- function(annovar_file, af_colname){

  
  #Remove all rows with LRG
  rm_lrg <- annovar_file %>%
    dplyr::filter(!grepl('LRG', Chr))
  
  #Select features needed for ML/EDA
  annovar_features <- rm_lrg %>%
    select(ID,
          Func.refGene,
           ExonicFunc.refGene, 
           af_colname,
           FATHMM_score,
           GenoCanyon_score,
           LRT_score,
           CADD_phred,
           `GERP++_RS`,
           `M-CAP_score`,
           MetaLR_score,
           MetaSVM_score, 
           MutationTaster_score,
           MutationAssessor_score,
           PROVEAN_score, 
           SIFT_score,
           SiPhy_29way_logOdds, 
           `fathmm-MKL_coding_score`,
           integrated_fitCons_score, 
           integrated_confidence_value, 
           phastCons100way_vertebrate,
           phastCons30way_mammalian,
           phyloP100way_vertebrate,
           phyloP30way_mammalian, 
           Otherinfo6) 
  
  
  #Replace "." with NA
  annovar_features[annovar_features == "."] <- NA
  annovar_features[annovar_features == ""] <- NA
  
  #Numeric features
  features_add <- annovar_features %>% 
    mutate_at(c(af_colname,
                'FATHMM_score',
                'GERP++_RS',
                'GenoCanyon_score',
                'LRT_score',
                "CADD_phred",
                'M-CAP_score',
                'MetaLR_score', 
                'MetaSVM_score', 
                'MutationTaster_score',
                'MutationAssessor_score',
                'PROVEAN_score', 
                'SIFT_score', 
                'SiPhy_29way_logOdds',
                'fathmm-MKL_coding_score',
                'integrated_fitCons_score', 
                'integrated_confidence_value',
                'phastCons100way_vertebrate',
                'phastCons30way_mammalian', 
                'phyloP100way_vertebrate',
                'phyloP30way_mammalian'),
              as.numeric) %>%
              rename(AF = af_colname)
  

}
rm_missing_rows <- function(df){
  
  nas <- df %>% 
    mutate(n_NAs = rowSums(is.na(.)))
  rm_missing <- subset(nas, n_NAs <= 18) %>%
    select(-n_NAs) %>%
    mutate(AF_NA = ifelse(is.na(AF), 1, 0)) %>%
    mutate(AF = replace_na(AF, 0)) 
  
}


# Preprocessing -----------------------------------------------------------

input_feat <- get_annotations(input, "${params.af}") %>%
  rm_missing_rows() %>%
  rename(`GERP.._RS` = `GERP++_RS`,
         `M.CAP_score` = `M-CAP_score`,
         `fathmm.MKL_coding_score` = `fathmm-MKL_coding_score`)


# Predict -----------------------------------------------------------------
if(nrow(input_feat) > 0) {
  preds <- add_predictions(input_feat,
                          mod,
                          type = "prob") %>%
    unnest(pred) %>%
    select(ID, `.pred_P_LP`)

  final_preds <- input %>%
    merge(preds,
          by = "ID") %>%
    arrange(desc(`.pred_P_LP`)) 

  write_tsv(final_preds,
    file = 'prioritised_file.txt')
} else {
  write.csv(input_feat,
            "prioritised_file.txt",
            row.names = FALSE)
}
  """
}


process SORT {
  publishDir params.outdir, mode:'copy'

  input:
  path merged_ml
  
  output:
  path 'prioritised_file.txt'
  
  script:
  """
  num_columns=\$(head -n 1 $merged_ml | tr '\t' '\n' | wc -l)
  { head -n 1 $merged_ml && tail -n +2 $merged_ml | sort -t '\t' -k\$num_columns,\$num_columns -rn; } > 'prioritised_file.txt'
  """
}



process TRIM {
  publishDir params.outdir, mode:'copy'
  
  input:
  path sorted_file
  
  output:
  path 'trimmed_prioritised_file.txt'
  
  script:
  """
  awk -F'\t' 'NR==1 || \$NF > 0.50' $sorted_file > 'trimmed_prioritised_file.txt'
  """
}


workflow{
    split_ch = SPLIT(params.vcf_file,
    params.indx_file) 
    filter_ch = FILTER(split_ch.flatten())
    annotate_ch = ANNOTATE(filter_ch).flatten()
    .filter( ~/(?:[^,]*\.txt)(?=,|$)/)
    ml_ch = ML(params.mod ,annotate_ch)
    .collectFile(name: 'merge.txt', 
    storeDir: params.outdir,
    keepHeader: true)
    sort_ch = SORT(ml_ch)
    trim_ch = TRIM(sort_ch)
}

workflow.onComplete = {
    println ( workflow.success ? "\nDone! Prioritised variants in $params.outdir/prioritised_file.txt\n" : "Oops .. something went wrong" )

}
