// main.nf
process JSON_TO_GTF {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*.gtf"), emit: finished_insert_overlap_transcript_gtf


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def json_to_gtf(json_path, output_gtf_path):
        with open(json_path, 'r') as json_file:
            gene_data = json.load(json_file)

        with open(output_gtf_path, 'w') as gtf_file:
            for gene_id in gene_data:
                for line in gene_data[gene_id]:
                    gtf_file.write(line + '\\n')

    # Generate output GTF file name based on JSON file
    output_gtf_name ='finished_insert_overlap_transcript.gtf'

    json_to_gtf('${json}', output_gtf_name)
    """
}

// workflow {
//     // Define the input file path
//     ch_json = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/updated_gene_id_data.json')

//     // Process the file
//     output_gtf = JSON_TO_GTF(ch_json.map { [ [:], it ] })
// }
