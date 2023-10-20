<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**mo-core/ifunmap** is a pipeline for Integrative FunMap network analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary
* FunMap generation
* Network analysis (degree, betweenness etc.)
* Dark gene enrichment analysis
* Network module identification ([ICE](http://ice.zhang-lab.org))
* Network module pathway enrichment analysis
* Network module visualization
* Module activity analysis (optional)
* Module activity prediction from mutation data (optional)


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Run QC analysis only.
   ```bash
   nextflow run mo-core/ifunmap -r main -latest --outdir <OUTDIR> --config_file /path/to/config_file.yml --data_file /path/to/data.tar.gz  -profile docker --qc_only
   ```

3. Run the FunMap network prediction and other network analysis steps:

   ```bash
   nextflow run mo-core/ifunmap -r main -latest --outdir <OUTDIR> --config_file /path/to/config_file.yml --data_file /path/to/data.tar.gz  -profile docker
   ```

  Here `config_file` and `data_file` are required. The default value for `outdir` is `results` under
  the current directory.

#### Instructions for Preparing configuration file for `funmap`

To configure a YAML specification file for your project, follow these steps:

*  **File Structure**: Structure your YAML file with key-value pairs. Use a colon (`:`) followed by a space to separate keys and values. You can include comments using the `#` symbol.

* **Define Task and Project Information**:
   - `task`: Specify the task type (`protein_func`).
   - `name`: Provide a unique name for your project.
   - `seed`: Set the seed value for randomization.

* **Results Directory and Filtering**:
   - `results_dir`: Define the directory where results will be stored.
   - `filter_noncoding_genes`: Set to `true` or `false` to enable or disable filtering of noncoding genes.

* **Feature Type and Job Settings**:
   - `feature_type`: Choose the feature type (`mr` for mutual rank or `cc` for Pearson correlation coefficient).
   - `n_jobs`: Specify the number of parallel jobs for tasks. (Defaults to all available cores)

* **Minimum Sample Count and Edge Parameters**:
   - `min_sample_count`: Set the minimum number of valid data points when computing correlations.
   - `start_edge_num`: Define the starting edge number.
   - `max_num_edges`: Specify the maximum edge number.
   - `step_size`: Set the step size for edge calculations.
   - `lr_cutoff`: Define the likelihood ratio cutoff to select valid edges in funmap.

* **Data Path and Data Files**:
   - `data_path`: Provide the path to the data archive.
   - `data_files`: List data files with the following information for each:
     - `name`: Name of the data file.
     - `type`: Type of the data file (e.g., 'protein' or 'rna').
     - `path`: Path to the data file.

* **RNA-Protein Pairs** (optional)
   - `rp_pairs`: Define RNA-protein pairs with the following details for each pair:
     - `name`: Name of the pair.
     - `rna`: Name of the corresponding RNA data file.
     - `protein`: Name of the corresponding protein data file.
