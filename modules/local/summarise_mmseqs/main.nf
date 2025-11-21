process SUMMARISE_MMSEQS {
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b2e8a2e5ca723da3a02451b91eea8219ac011a303d01a4f2b1f15ea58fb79fb5/data'
        : 'community.wave.seqera.io/library/polars:1.35.1--f9b832ec070b095b'}"

    input:
    tuple val(meta), path(mmseqs_tsv)

    output:
    path ("*.csv"), emit: summary_csv
    path ("taxlineage/${prefix}_taxlineage.csv"), emit: taxlineage
    tuple val(meta), path("reads_*.tsv"), emit: abundance_picrust
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: implement this better
    metadata = "sample=${meta.id}"
    if (meta.containsKey('group') && meta.group != 'false') {
        metadata += ",group=${meta.group}"
    }
    """
    summarise_mmseqs.py \\
        --mmseqs-tsv ${mmseqs_tsv} \\
        --min-aln ${params.mmseqs2_min_aln} \\
        --min-identity ${params.mmseqs2_min_identity} \\
        --metadata ${metadata}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir taxlineage
    touch taxlineage/${prefix}_taxlineage.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        polars: \$(python -c "import polars; print(polars.__version__)")
    END_VERSIONS
    """
}
