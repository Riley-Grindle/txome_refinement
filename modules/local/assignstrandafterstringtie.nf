process ASSIGN_STRAND_AFTER_STRINGTIE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::sed=4.7 conda-forge::grep=3.11 conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"


    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("strand_after_StringTie.gtf") , emit: processed_gtf
    tuple val(meta), path("delete_strand_after_StringTie.gtf"), emit: deleted_gtf
    tuple val(meta), path("only_modified_strand_part_after_StringTie.gtf"), emit: modified_gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash

    touch delete_strand_after_StringTie.gtf
    touch only_modified_strand_part_after_StringTie.gtf

    awk -F'\t' '{
        if (\$7 == ".") {
            n = split(\$9, attrs, ";"); # Split attributes
            for (i = 1; i <= n; i++) {
                if (match(attrs[i], /cov "([^"]+)"/)) {
                    cov_value = substr(attrs[i], RSTART + 5, RLENGTH - 6); # Extract cov value
                    cov_value += 0; # Convert to number
                    if (cov_value < 10) {
                        print \$0 >> "delete_strand_after_StringTie.gtf"; # Record deleted lines
                        next; # Skip lines with coverage less than 10
                    } else {
                        \$7 = "+"; # Set strand to '+'
                        print \$0 >> "only_modified_strand_part_after_StringTie.gtf"; # Record lines modified to '+'
                    }
                    break; # Exit loop after finding cov
                }
            }
        }
        print \$0;
    }' OFS='\t' "${gtf}" > "strand_after_StringTie.gtf"
    """
}
