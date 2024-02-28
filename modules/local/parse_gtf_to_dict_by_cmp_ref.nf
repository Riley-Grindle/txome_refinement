process PARSE_GTF_TO_DICT {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.json"), emit: parsed_gtf_json_by_cmp_ref

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def parse_gtf_to_dict_by_cmp_ref(file_path):
        cmp_ref_data = {}
        with open(file_path, 'r') as file:
            current_cmp_ref = None
            for line in file:
                if line.startswith("#"):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                if fields[2] == "transcript" and 'cmp_ref' in fields[8]:
                    cmp_ref_field = [f for f in fields[8].split(';') if 'cmp_ref' in f]
                    if cmp_ref_field:
                        current_cmp_ref = cmp_ref_field[0].split('"')[1]
                        if current_cmp_ref not in cmp_ref_data:
                            cmp_ref_data[current_cmp_ref] = []

                if current_cmp_ref:
                    cmp_ref_data[current_cmp_ref].append(line.strip())

        return cmp_ref_data

    # Generate output JSON file name based on GTF file
    output_json_name = os.path.basename('${gtf}').replace('.gtf', '_cmp_ref_data.json')

    cmp_ref_dict = parse_gtf_to_dict_by_cmp_ref('${gtf}')

    with open(output_json_name, 'w') as json_file:
        json.dump(cmp_ref_dict, json_file, indent=4)
    """
}

// workflow {
//     // Define the input file path
//     ch_gtf = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/readyInsert.gtf')

//     // Process the file
//     parsed_gtf_json = PARSE_GTF_TO_DICT(ch_gtf.map { [ [:], it ] })
// }
