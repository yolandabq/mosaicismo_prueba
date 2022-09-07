process VARSCAN_MPILEUP2SNP {
    
    tag "$meta.id"
    label 'process_low'
    conda (params.enable_conda ? "bioconda::varscan=2.4.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.4--hdfd78af_1':
        'quay.io/biocontainers/varscan:2.4.4--hdfd78af_1' }"
        
    input:
    tuple val(meta), path(mpileup)
    
    output:									// It defines the output of the process (i.e. files) and send to a new channel
    path "*.vcf", emit: vcf_mpileup
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
        varscan mpileup2snp ${mpileup} --min-var-freq 0.2 --p-value 1 --output-vcf 1 > ${prefix}_varscan.vcf
        cat <<-END_VERSIONS > versions.yml
    """
}

