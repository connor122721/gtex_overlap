#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import the required modules
include { runColoc } from './modules/coloc.nf'
include { analysisColoc } from './modules/coloc.nf'
include { execute_replication } from './modules/replication_eQTLs.nf'

// GTEx data: extract chromosome from file name
gtex = Channel.fromPath("${params.gtex_dir}/Heart_Left_Ventricle.v10.allpairs.*.parquet")
gtex_with_chr = gtex.map { file ->
    def fileName = file.getFileName().toString()
    def m = fileName =~ /chr(\w+)\.parquet$/
    if(m) { [ "chr" + m[0][1], file.toAbsolutePath() ] }
    else { error("Filename ${fileName} does not match expected format") }
}
//gtex_with_chr.view()

// GWAS data
gwas = Channel.fromPath("${params.gwas_dir}/processed_jurgens24*.rds")
gwas_with_chr = gwas.map { file ->
    def fileName = file.getFileName().toString()
    def m = fileName =~ /chr(\w+)\.rds$/
    if(m) { [ "chr" + m[0][1], file.toAbsolutePath() ] }
    else { error("Filename ${fileName} does not match expected format") }
}
//gwas_with_chr.view()

// Combine by chromosome
gtex_with_chr
    .combine(gwas_with_chr, by: 0)
    .map { [ it[1], it[0], it[2] ] }
    .set { gtex_gwas }
//gtex_gwas.view()

// Run Workflow!
workflow {

    // Replication analyses
    // execute_replication(params.gtex_dir, params.chef_file, params.output_dir)

    // Parameters for coloc analyses
    N_jurgens = Channel.of(955733, 516).toList()

    // 1) Run coloc
    colocJurgens = runColoc(gtex_gwas
                            .combine(N_jurgens))

    // Get candidate eGenes that are colocalized and prep for LD analysis
    wd1 = colocJurgens.outDir.unique().collect() 
    ColocGenes = analysisColoc(wd1)
    
}
