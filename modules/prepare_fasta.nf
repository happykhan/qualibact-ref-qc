process PREPARE_FASTA {
    tag { meta.sample_id }
    label 'process_mini'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.sample_id}.fasta")

    script:
    def out = "${meta.sample_id}.fasta"

    """
    set -euo pipefail
    echo "[PREPARE_FASTA] meta.sample_id=${meta.sample_id} fasta='${fasta}' out='${out}'"

    # Try GNU readlink -f first, otherwise use a small python fallback.
    # Use a single command (no heredoc) so substitution works reliably.
    if readlink -f "${fasta}" >/dev/null 2>&1; then
        inpath=\$(readlink -f "${fasta}")
    else
        inpath=\$(python -c 'import os,sys; print(os.path.realpath(sys.argv[1]))' "${fasta}")
    fi

    # out may not exist yet — make a canonical path for comparison
    if readlink -f "${out}" >/dev/null 2>&1; then
        outpath=\$(readlink -f "${out}")
    else
        outpath="\$PWD/${out}"
    fi

    echo "[PREPARE_FASTA] inpath=\$inpath"
    echo "[PREPARE_FASTA] outpath=\$outpath"

    # If input is gzipped -> decompress to tmp then move into place
    if [[ "${fasta}" == *.gz ]]; then
        tmp="\${out}.tmp"
        echo "[PREPARE_FASTA] detected gz input; writing to \$tmp"
        gzip -cd "${fasta}" > "\$tmp"
        mv -f "\$tmp" "${out}"
        echo "[PREPARE_FASTA] decompressed -> ${out}"
        exit 0
    fi

    # If input and output are the SAME file -> nothing to do
    if [[ "\$inpath" == "\$outpath" ]]; then
        echo "[PREPARE_FASTA] input and output are identical; nothing to do"
        # ensure Nextflow will see the path
        if [[ ! -e "${out}" ]]; then
            ln -s "\$inpath" "${out}"
            echo "[PREPARE_FASTA] created symlink ${out} -> \$inpath"
        fi
        exit 0
    fi

    # Otherwise copy to tmp then atomically move into place
    tmp="\${out}.tmp"
    echo "[PREPARE_FASTA] copying ${fasta} -> \$tmp"
    cp -p "${fasta}" "\$tmp"
    mv -f "\$tmp" "${out}"
    echo "[PREPARE_FASTA] wrote ${out}"
    """
}
