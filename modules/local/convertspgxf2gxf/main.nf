
process AGAT_CONVERTSPGXF2GXF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.0.0--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    path("*.agat.gtf")  , emit: output_gtf, optional: true
    path("*.final.gtf")  , emit: final_gtf, optional: true
    path("*.log")       , emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: 'agat'
    """
    agat config --expose

    sed -i.bak 's/output_format: *.*/output_format: GTF/' agat_config.yaml
    rm *.bak

    agat_convert_sp_gxf2gxf.pl \\
        -g $gtf \\
        -o ${prefix}.${args2}.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gxf2gxf.pl --help | sed '4!d; s/.*v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.agat.gff
    touch ${gff}.agat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_convert_sp_gxf2gxf.pl --help | sed '4!d; s/.*v//')
    END_VERSIONS
    """
}

