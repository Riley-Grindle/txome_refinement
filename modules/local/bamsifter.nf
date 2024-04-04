process BAMSIFTER {
    tag "$meta.id"
    label 'process_bamsifter'


    container "docker.io/trinityctat/ctat_vif"

    input:
    tuple val(meta), path(alignment_bam)

    output:
    tuple val(meta), path ("*.bam")                     , emit: normalized_bam
    //path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when
    alignment_bam.size() > 0

    script:

    def args = task.ext.args ?: ''

    def normalize_max_read_cov = 30  // better for polymorphic transcriptomes

    //def normalize_max_read_cov = 200  // better for polymorphic transcriptomes
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /usr/local/src/CTAT-VirusIntegrationFinder/util/bamsifter/bamsifter \\
            -c ${normalize_max_read_cov} \\
            -o ${prefix}.${args}sifted.bam \\
            ${alignment_bam}
    """

}
