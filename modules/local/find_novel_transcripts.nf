// main.nf
process FIND_NOVEL_TRANSCRIPTS {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta), path(reference_gtf)

    // output:
    // tuple val(meta), path("*.gtf"), emit: select_transcript_append_based_on_gffcompare_class_code_to_generate_final_annotation_gtf
    output:
    tuple val(meta), path("*.gtf"), emit: final_annotation_gtf

    when:
    task.ext.when == null || task.ext.when



    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    # Convert the JSON string to a Python list
    novel_class_codes = '${params.novelClassCodes}'

    def find_novel_transcripts(gtf, reference_gtf_path, novel_class_codes):


        # Generate output file name based on combined GTF file
        output_gtf_name ='final_annotation.gtf'

        # Read the combined GTF file
        with open(gtf, 'r') as file:
            combined_gtf_lines = file.readlines()

        novel_transcript_lines = []
        target_class_code = False
        for line in combined_gtf_lines:
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            # Check for transcript lines with novel class codes
            if fields[2] == "transcript" and 'class_code "' in fields[8]:
                class_code = fields[8].split('class_code "')[1][0]
                target_class_code = class_code in novel_class_codes
                if target_class_code:
                    novel_transcript_lines.append(line)

            elif target_class_code:
                novel_transcript_lines.append(line)

        # Read the reference GTF file and append the novel transcript lines
        with open(reference_gtf_path, 'r') as ref_file:
            reference_gtf_lines = ref_file.readlines()

        # Write the combined reference and novel transcript lines to the output GTF file
        with open(output_gtf_name, 'w') as output_file:
            output_file.writelines(reference_gtf_lines + novel_transcript_lines)

        return output_gtf_name

    output_gtf_path = find_novel_transcripts('${gtf}', '${reference_gtf}',novel_class_codes)
    """
}

// workflow {
//     // Define the input file paths
//     ch_gtf = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/null.combined.gtf')
//     ch_ref_gtf = Channel.fromPath('/Users/dxu/Documents/axolotlOmics/AmexT_v47-AmexG_v6.0-DD.gtf')

//     // Process the file
//     output_gtf_path = FIND_NOVEL_TRANSCRIPTS(ch_gtf.map { [ [:], it ] },ch_ref_gtf.map { [ [:], it ]} )
// }
