#!/usr/bin/env nextflow

params.gtex_dir = "/standard/vol185/cphg_Manichaikul/users/csm6hg/gtex_overlap/data/gtex_full_ciseqtl"
params.chef_file = "/standard/vol185/cphg_Manichaikul/users/csm6hg/data/candidate_maxRNAPC11_eQTL_sigPerm0.05.rds"
params.output_dir = "/standard/vol185/cphg_Manichaikul/users/csm6hg/gtex_overlap/results"

process execute_replication {
    publishDir params.output_dir, mode: 'copy'
    
    input:
    val gtex_dir
    val chef_file
    val output_dir
    
    output:
    path "overlapGTEX_eqtl_results.csv"
    path "overlap_plot.png"
    
    script:
    """
    python ${projectDir}/bin/gtex_replication.py \
        --gtex_dir ${gtex_dir} \
        --chef_file ${chef_file} \
        --output_dir ${output_dir}
    """
}

workflow {
    execute_replication(params.gtex_dir, params.chef_file, params.output_dir)
}
