// main.nf
process INSERT_BY_START_POSITION {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(cmp_ref_data)
    tuple val(meta), path(reference_gtf)

    output:
    tuple val(meta), path("*.json"), emit: finished_insert_overlap_transcript_json

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    def separate_transcripts_with_same_cmp_ref_geneID(gtf_list):
        processed_dict = {}
        key_counter = 1
        for line in gtf_list:
            fields = line.split('\t')
            if fields[2] in ['transcript', 'exon']:
                start_position = fields[3] if fields[2] == 'transcript' else start_position
                processed_dict.setdefault(start_position, []).append(line)
            else:
                unique_key = f'unique_key_{key_counter}'
                processed_dict[unique_key] = [line]
                key_counter += 1
        return processed_dict

    def find_last_transcript_key(reference_value):
        last_transcript_key = None
        for key in reference_value:
            if key.startswith("unique_key_"):
                continue
            last_transcript_key = key
        return last_transcript_key

    # Load JSON files
    with open('${cmp_ref_data}', 'r') as file:
        cmp_ref_data = json.load(file)

    with open('${reference_gtf}', 'r') as file:
        reference_data = json.load(file)

    # Process each cmp_ref_data
    for cmp_ref, lines in cmp_ref_data.items():
        if cmp_ref in reference_data:
            cmp_ref_value = separate_transcripts_with_same_cmp_ref_geneID(lines)
            reference_value = separate_transcripts_with_same_cmp_ref_geneID(reference_data[cmp_ref])

            for cmp_ref_start_position, cmp_ref_lines in cmp_ref_value.items():
                inserted = False
                for reference_start_position, reference_lines in reference_value.items():
                    if reference_start_position.isdigit() and int(cmp_ref_start_position) <= int(reference_start_position):
                        reference_value[reference_start_position] = cmp_ref_lines + reference_lines
                        inserted = True
                        break

                if not inserted:
                    last_transcript_key = find_last_transcript_key(reference_value)
                    if last_transcript_key is not None:
                        reference_value[last_transcript_key].extend(cmp_ref_lines)
                    else:
                        unique_key = f'unique_end_key_{cmp_ref_start_position}'
                        reference_value[unique_key] = cmp_ref_lines

            reference_data[cmp_ref] = [line for sublist in reference_value.values() for line in sublist]

    # Generate output JSON file name
    output_json_name = 'finished_insert_overlap_transcript.json'
    with open(output_json_name, 'w') as json_file:
        json.dump(reference_data, json_file, indent=4)
    """
}

// workflow {
//     // Define the input file paths
//     ch_cmp_ref_data = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/cmp_ref_data_after_replace_key_trackingFile.json')
//     ch_reference_gtf = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/gene_id_data.json')

//     // Process the file
//     updated_reference_data_json = INSERT_BY_START_POSITION(ch_cmp_ref_data.map { [ [:], it ] }, ch_reference_gtf.map { [ [:], it ] })
// }
