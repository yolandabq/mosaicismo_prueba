process VARDICTJAVA {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::vardict-java=1.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0':
        'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(input) //path(bam), path(bai), path(bed)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.8.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    vardict-java \\
        $args \\
        -c 1 -S 2 -E 3 \\
        -b ${input[0]} \\
        -th $task.cpus \\
        -N $prefix \\
        -G $fasta \\
        ${input[2]} \\
        | teststrandbias.R \\
        | var2vcf_valid.pl \\
            $args2 \\
            -N $prefix > ${prefix}_vardictjava.vcf #\\
        #| gzip -c > ${prefix}_vardictjava.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: $VERSION
        var2vcf_valid.pl: \$(echo \$(var2vcf_valid.pl -h | sed -n 2p | awk '{ print \$2 }'))
    END_VERSIONS
    """
}
