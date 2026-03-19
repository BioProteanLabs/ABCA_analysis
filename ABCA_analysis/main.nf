#!/usr/bin/env nextflow

// Define parameters for scripts and output directories
params.data_access    = '01_data_access.ipynb'
params.add_meta       = '02_add_meta_to_expr.ipynb'
params.labels_plots   = '03_labels_plots.ipynb'
params.me_explore     = '04_ME_explore_plots.ipynb'
params.partial_mods   = '05_ME_partial_modules.ipynb'
params.me_tau         = '06_ME_tau.ipynb'
params.me_tau_table   = '07_ME_tau_table.py'
params.cell_typing    = '08_cell_typing.ipynb'

params.output_dir     = '/data/scRNA/ABCA/'
params.data_dir       = 'data'

// Ensure output directory exists
process setup {
    output:
        path params.output_dir

    script:
    """
    mkdir -p ${params.output_dir}
    """
}

// Process 1: Data access
process dataAccess {
    input:
        path setup_done from setup
    output:
        path "${params.output_dir}/data_out.h5ad" into data_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.data_access} ${params.output_dir}/01_data_access_output.ipynb
    mv path/to/generated_data.h5ad ${params.output_dir}/data_out.h5ad
    """
}

// Process 2: Add metadata to expression data
process addMeta {
    input:
        path data from data_ch
    output:
        path "${params.output_dir}/data_with_meta.h5ad" into meta_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.add_meta} ${params.output_dir}/02_add_meta_output.ipynb -p input_data ${data}
    mv path/to/modified_data.h5ad ${params.output_dir}/data_with_meta.h5ad
    """
}

// Process 3: Labels and plots
process labelsPlots {
    input:
        path meta from meta_ch
    output:
        path "${params.output_dir}/labels_plots_output.png" into labels_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.labels_plots} ${params.output_dir}/03_labels_plots_output.ipynb -p input_data ${meta}
    mv path/to/labels_plots_output.png ${params.output_dir}/labels_plots_output.png
    """
}

// Process 4: Module Eigengene Exploration
process exploreME {
    input:
        path meta from meta_ch
    output:
        path "${params.output_dir}/me_explore_output.pdf" into me_explore_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.me_explore} ${params.output_dir}/04_ME_explore_plots_output.ipynb -p input_data ${meta}
    mv path/to/me_explore_output.pdf ${params.output_dir}/me_explore_output.pdf
    """
}

// Process 5: Partial Module Analysis
process partialModules {
    input:
        path meta from meta_ch
    output:
        path "${params.output_dir}/partial_modules_output.txt" into partial_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.partial_mods} ${params.output_dir}/05_ME_partial_modules_output.ipynb -p input_data ${meta}
    mv path/to/partial_modules_output.txt ${params.output_dir}/partial_modules_output.txt
    """
}

// Process 6: Compute ME Tau using notebook
process meTau {
    input:
        path meta from meta_ch
    output:
        path "${params.output_dir}/tau_output.txt" into tau_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.me_tau} ${params.output_dir}/06_ME_tau_output.ipynb -p input_data ${meta}
    mv path/to/tau_output.txt ${params.output_dir}/tau_output.txt
    """
}

// Process 7: Compute Tau Table (Python Script)
process meTauTable {
    input:
        path meta from meta_ch
    output:
        path "${params.output_dir}/tau_table.csv" into tau_table_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    python ${params.me_tau_table}
    mv /data/scRNA/ABCA/AIBS/AWS/expression_matrices/WMB-10Xv3/20230630/outputs/WMB-10Xv3-Isocortex-1-raw-sc-wgcna-ME-sc-mergedTau.csv ${params.output_dir}/tau_table.csv
    """
}

// Process 8: Cell Typing Analysis
process cellTyping {
    input:
        path meta from meta_ch
    output:
        path "${params.output_dir}/cell_typing_output.txt" into cell_typing_ch

    conda 'envs/abca-env.yaml'
    script:
    """
    papermill ${params.cell_typing} ${params.output_dir}/08_cell_typing_output.ipynb -p input_data ${meta}
    mv path/to/cell_typing_output.txt ${params.output_dir}/cell_typing_output.txt
    """
}

// Define the workflow connections.
workflow {
    setup_done = setup()
    data = dataAccess(setup_done)
    meta = addMeta(data)
    labelsPlots(meta)
    exploreME(meta)
    partialModules(meta)
    meTau(meta)
    meTauTable(meta)
    cellTyping(meta)
}
