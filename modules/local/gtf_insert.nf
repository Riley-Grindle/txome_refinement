
process GTF_INSERT {

    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/gtf_insert' }"

    
    input:
    tuple val(meta), path(combined_gtf)
    tuple val(meta), path(tracking_file)
    tuple val(meta), path(reference_gtf)

    output:
    path("final_annotation.gtf"), emit: final_gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    find_overlaps.py $combined_gtf
    parse_gtf_to_dict_cmp_ref.py find_overlaps.gtf
    extract_key_value.py $tracking_file
    check_unique_keys.py key_value_genes.json
    update_keys.py key_value_genes.json gtf_to_dict_gffcmp.json
    parse_gtf_2_dict_geneid.py $reference_gtf
    insert_by_start.py updated_gffcmp.json parse_reference_gtf.json
    json_2_gtf.py overlap_inserted.json
    find_novel_tscripts.py $combined_gtf overlap_inserted.gtf $tracking_file
    """
}

