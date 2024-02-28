process FIND_OVERLAP_TRANSCRIPT {

    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gtf"), emit: select_transcript_insert_based_on_gffcompare_class_code

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import os
    import json

    # Assuming overlapClassCodes is passed as '["k", "m", "n", "j", "e"]'
    overlap_class_codes = '${params.overlapClassCodes}'


    def find_overlap_transcript(gtf, overlap_class_codes):
        # Generate output file name
        output_gtf_name = os.path.basename(gtf).replace('.gtf', '_ready_insert_transcript.gtf')

        # Open the input GTF file
        with open(gtf, 'r') as file:
            combined_gtf_lines = file.readlines()

        overlap_transcript_lines = []
        for line in combined_gtf_lines:
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] == "transcript" and 'class_code "' in fields[8]:
                class_code = fields[8].split('class_code "')[1][0]
                if class_code in overlap_class_codes:
                    overlap_transcript_lines.append(line)
                    # Add all subsequent exon lines
                    index = combined_gtf_lines.index(line) + 1
                    while index < len(combined_gtf_lines) and 'exon' in combined_gtf_lines[index].split('\t')[2]:
                        overlap_transcript_lines.append(combined_gtf_lines[index])
                        index += 1

        # Write the overlap transcript lines to the output GTF file
        with open(output_gtf_name, 'w') as output_file:
            output_file.writelines(overlap_transcript_lines)

        return output_gtf_name

    # Call the function and pass the list of overlap class codes
    output_gtf_path = find_overlap_transcript('${gtf}', overlap_class_codes)
    """
}



// workflow {
// // Define the input file path
// ch_gtf = Channel.fromPath('/Users/dxu/Documents/compareJoelGTFwithMyGTF/updateGeneID/outputFromGffcompare/null.combined.gtf')

// // Process the file
// output_gtf_path = FIND_OVERLAP_TRANSCRIPT(ch_gtf.map { [ [:], it ] })
// }
