// main.nf
process PARSE_GTF_TO_DICT_BY_GENEID {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.json"), emit: parsed_gtf_json_by_geneID

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def parse_gtf_to_dict_by_geneID(file_path):
        gene_data = {}
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue  # Skip comment lines
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                gene_id_field = [f for f in fields[8].split(';') if 'gene_id' in f]
                if gene_id_field:
                    gene_id = gene_id_field[0].split('"')[1]

                    if gene_id not in gene_data:
                        gene_data[gene_id] = []

                    gene_data[gene_id].append(line.strip())

        return gene_data

    # Generate output JSON file name based on GTF file
    output_json_name = os.path.basename('${gtf}').replace('.gtf', '_gene_id_data.json')

    gene_id_dict = parse_gtf_to_dict_by_geneID('${gtf}')

    with open(output_json_name, 'w') as json_file:
        json.dump(gene_id_dict, json_file, indent=4)
    """
}

// workflow {
//     // Define the input file path
//     ch_gtf = Channel.fromPath('/Users/dxu/Documents/axolotlOmics/AmexT_v47-AmexG_v6.0-DD.gtf')

//     // Process the file
//     parsed_gtf_json = PARSE_GTF_TO_DICT_BY_GENEID(ch_gtf.map { [ [:], it ] })
// }
