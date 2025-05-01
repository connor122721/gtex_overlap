// Replication analyses
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