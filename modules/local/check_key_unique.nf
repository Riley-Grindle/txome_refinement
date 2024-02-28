// main.nf
process CHECK_UNIQUE_KEYS {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(json_file)

    output:
    tuple val(meta), path("results.txt"), emit: unique_keys_result

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def check_unique_keys(json_file):
        try:
            with open(json_file, 'r') as file:
                data = json.load(file)

            keys = list(data.keys())
            if len(keys) == len(set(keys)):
                return True
            else:
                return False
        except Exception as e:
            return str(e)

    result = check_unique_keys('${json_file}')
    result_text = f"Are all keys unique? {result}"

    # Save the result to a text file
    output_text_name = 'results.txt'
    with open(output_text_name, 'w') as output_file:
        output_file.write(result_text)
    """

}

// workflow {
//     // Define the input file path
//     ch_json_file = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/keyvalue.json')

//     // Process the file
//     unique_keys_result = CHECK_UNIQUE_KEYS(ch_json_file.map { [ [:], it ] })
// }
