#!/usr/bin/env nextflow

//define parameters
 params.vcf_file = "$projectDir/vcf_test/madcap_test2.vcf.gz"
 params.indx_file = "$projectDir/vcf_test/madcap_test2.vcf.gz.csi"
 params.annovar_db = "/home/ckitching/lustre/anno/humandb"
 params.mod = "$projectDir/rf_mod-04-23.rds"
 params.rscript = "$projectDir/predictions.R"
 params.split_script = "$projectDir/split.sh"
 params.build = "hg19"
params.outdir = "results"

/*conditional gnomad params: > gnomad 2 NA to hg19
    default version for hg38 is gnomad 4
*/
if ("${params.build}" == "hg19" ) {
    println "Using hg19."
    params.gnomad = "gnomad211_genome"
}
else if ("${params.build}" == "hg38" ) {
    params.gnomad = "gnomad40_genome"
}


 log.info """\
    V P   P I P E L I N E
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
    bash ${params.split_script} $vcf $indx chr
    
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
  publishDir params.outdir, mode:'copy'

  input:
  path anno_file
  
  output:
  path 'prioritised_file.txt'
  
  script:
  """
  Rscript ${params.rscript} ${params.mod} $anno_file prioritised_file.txt

  """
}



workflow{
    split_ch = SPLIT(params.vcf_file,
    params.indx_file) 
    filter_ch = FILTER(split_ch.flatten())
    annotate_ch = ANNOTATE(filter_ch).flatten()
    .filter( ~/(?:[^,]*\.txt)(?=,|$)/ )
    .collectFile(name: 'merge.txt', 
    storeDir: params.outdir,
    keepHeader: true)
    annotate_ch.view()
    ml_ch = ML(annotate_ch)
}

workflow.onComplete = {
    println ( workflow.success ? "\nDone! Prioritised variants in $params.outdir/prioritised_file.txt\n" : "Oops .. something went wrong" )

}
