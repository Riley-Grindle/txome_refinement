process Staging {

    conda "bioconda::trinity=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.13.2--h00214ad_1':
        'biocontainers/trinity:2.13.2--h00214ad_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    path("samples_file.txt"), emit: samplesFileChannel

    script:
    def read1 = reads[0].toString()
    def read2 = reads.size() > 1 ? reads[1].toString() : ""
    def sample_entry = [meta.id, "${meta.id}_${meta.strandedness}", read1, read2].join('\t')

    """
    echo "$sample_entry" >> samples_file.txt

    """
}
// process Staging {
//     input:
//     tuple val(meta), path(reads, stageAs: "input*/*")

//     output:
//     path("samples_single_end.txt"), emit: single_end_samples
//     path("samples_double_end.txt"), emit: double_end_samples

//     script:
//     def read1 = reads[0].toString()
//     def read2 = reads.size() > 1 ? reads[1].toString() : ""
//     def sample_entry = [meta.id, "${meta.id}_${meta.strandedness}", read1, read2].join('\t')

//     if(meta.single_end == true) {
//         """
//         echo "$sample_entry" >> samples_single_end.txt
//         """
//     } else {
//         """
//         echo "$sample_entry" >> samples_double_end.txt
//         """
//     }
// }
