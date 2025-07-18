# `moducomp`: metabolic module completeness of genomes and metabolic complementarity in microbiomes

`moducomp` is a bioinformatics pipeline designed to identify and analyze metabolic module completeness and complementarity in microbial communities. It processes genomic data (protein sequences in FAA format) to map KEGG Orthology (KO) terms to KEGG modules, calculates module completeness for individual genomes and combinations of genomes within a microbiome, specifically reports potential complementary relationships where a metabolic module is 100% complete in a combination of N genomes but not in any individual member or smaller subset.

## Features

- Annotation of protein sequences using [`eggNOG-mapper`](https://github.com/eggnogdb/eggnog-mapper) to obtain KO terms.
- Mapping of KOs to KEGG metabolic modules based on [`kegg-pathways-completeness-tool`](https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool) to obtain metabolic module completeness.
- **Parallel processing support** for faster KPCT (KEGG Pathways Completeness Tool) analysis with automatic chunking and checkpointing.
- Reporting of module completeness for individual genomes.
- Calculation of module completeness for N-member genome combinations.
- Generation of complementarity reports highlighting modules completed through genome partnerships.
- Tracks and reports the actual proteins that are responsible for the completion of the module in the combination of N genomes.

## Installation

`moducomp` uses [Pixi](https://pixi.sh/) for managing dependencies and environment setup. Pixi ensures that you have a consistent and reproducible environment.

1.  **Install Pixi**:
If you don't have Pixi installed, follow the instructions on the [official Pixi website](https://pixi.sh/latest/#installation).

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

2.  **Clone the Repository** (if you haven't already):
```bash
git clone https://github.com/NeLLi-team/moducomp.git
cd moducomp
```

3.  **Install Dependencies using Pixi**:
Navigate to the project's root directory (where `pixi.toml` is located) and run:
```bash
pixi install
```
This command will read the `pixi.toml` file, resolve the dependencies (including Python, eggNOG-mapper, kegg-pathways-completeness-tool, pandas, etc.), and set up a local environment in a `.pixi` subdirectory.

*Note on eggNOG-mapper data*: download the EggNOG-mapper DB that is required for functional annotation. This can be a large download and may take time. Ensure you have sufficient disk space and an internet connection:

```bash
pixi shell
export EGGNOG_DATA_DIR="/path/to/datadir/to/store/eggnog-db"
download_eggnog_data.py
```

Alternatively, you can configure eggNOG-mapper to use a pre-downloaded database.
```bash
export EGGNOG_DATA_DIR="/path/to/datadir/to/store/eggnog-db"
```

## Usage

`moducomp` provides two main commands: `pipeline` and `analyze-ko-matrix`. You can run these commands using Pixi tasks defined in `pixi.toml` or directly within the Pixi environment.

### Performance and Parallel Processing

`moducomp` includes **parallel processing capabilities** for the KPCT (KEGG Pathways Completeness Tool) analysis, which can significantly improve performance for large datasets:

- **Automatic chunking**: Input files are automatically split into chunks for parallel processing
- **Checkpointing**: Resume capability if processing is interrupted - already processed chunks are automatically detected and skipped
- **Fallback mechanism**: If parallel processing fails, the system automatically falls back to sequential processing
- **Configurable CPU usage**: Use the `--ncpus` parameter to control how many CPU cores to use

**CPU Configuration**:
- The `--ncpus` parameter controls the number of CPU cores used for both eggNOG-mapper annotation and KPCT analysis
- For KPCT parallel processing, the system creates the same number of chunks as CPU cores specified
- Example: `--ncpus 8` will use 8 cores and create 8 chunks for optimal parallel processing

### Important Note

`moducomp` is specifically designed for large scale analysis of microbiomes with hundreds of members, and works on Linux systems with at least **64GB of RAM**. Nevertheless, it can be run on **smaller systems with less RAM, using the flag `--lowmem` when running the `pipeline` command**.

To activate the Pixi environment shell:
```bash
pixi shell
```
Once in the shell, you can run `./moducomp.py --help` for a full list of commands and options.

### Testing the Pipeline (Example)

To test the pipeline, you'll need a directory with some sample genome FAA files.

1.  **Prepare FAA files**: Create a directory (e.g., `test_genomes/`) and place your sample `.faa` files in it.
Ensure FAA headers are simple (e.g., `>protein_id`) or use the `--adapt-headers` option if they follow a `>genomeName|proteinId` format that needs processing into `>genomeName_prefix|protein_id_counter` for `emapper.py`.

2.  **Run the pipeline command**:

If you have enough RAM (>64GB), you can run the pipeline with the following command:

```bash
pixi shell
./moducomp.py pipeline ./test_genomes ./output_test_pipeline --ncpus 16 --calculate-complementarity 3
```

⚠️⚠️️⚠️ Note: If don't have enough RAM (>64GB), you can use the `--lowmem` flag to run the pipeline with less memory but it **will be slower**.

```bash
pixi shell
./moducomp.py pipeline ./test_genomes ./output_test_pipeline_lowmem --ncpus 16 --calculate-complementarity 3 --lowmem
```

### Running with your samples

#### `pipeline` command

Use the `pipeline` command to process a directory of genome FAA files from scratch.

```bash
pixi shell

./moducomp.py pipeline \
    /path/to/your/faa_files \
    /path/to/your/output_directory \
    --ncpus <number_of_cpus_to_use> \
    --calculate-complementarity <N>  # 0 to disable, 2 for 2-member, 3 for 3-member complementarity.
    # Optional flags:
    # --lowmem                    # Optional: Use this if you have less than 64GB of RAM
    # --adapt-headers             # If your FASTA headers need modification
    # --del-tmp                   # To delete temporary files
```

#### `analyze-ko-matrix` command

Use the `analyze-ko-matrix` command if you already have a KO matrix file (CSV format, where rows are genomes/combinations, columns are KOs, and values are KO counts).

The KO matrix file should have a `taxon_oid` column for genome identifiers, and subsequent columns for each KO (e.g., `K00001`, `K00002`) with integer counts.

```bash
pixi shell

./moducomp.py analyze-ko-matrix \
    /path/to/your/kos_matrix.csv \
    /path/to/your/output_directory \
    --ncpus <number_of_cpus_to_use> \
    --calculate-complementarity <N>  # 0 to disable, 2 for 2-member, 3 for 3-member complementarity.

    # Optional flags:
    # --del-tmp false
```

### Parallel Processing Features

`moducomp` includes advanced parallel processing capabilities for improved performance:

#### KPCT Parallel Processing

When using the `--ncpus` parameter with a value greater than 1, `moducomp` automatically enables parallel processing for the KPCT (KEGG Pathways Completeness Tool) analysis:

- **Automatic chunking**: Input files are split into `ncpus` chunks for optimal load balancing
- **Concurrent processing**: Multiple chunks are processed simultaneously using `multiprocessing.ProcessPoolExecutor`
- **Resume capability**: If processing is interrupted, completed chunks are automatically detected and skipped on restart
- **Automatic fallback**: If parallel processing fails, the system seamlessly falls back to sequential processing

#### Performance Tips

- **CPU cores**: Start with `--ncpus 8` for moderate datasets, increase to `--ncpus 16` or higher for large datasets
- **Memory considerations**: Each parallel worker requires memory; reduce `--ncpus` if you encounter memory issues
- **Large datasets**: For datasets with hundreds of genomes, parallel processing can reduce KPCT analysis time by 50-80%

#### Example with Parallel Processing

```bash
# For large datasets with sufficient resources
./moducomp.py pipeline ./large_genome_collection ./output_large --ncpus 32 --calculate-complementarity 3

# For moderate datasets
./moducomp.py analyze-ko-matrix ./ko_matrix.csv ./output_moderate --ncpus 16 --calculate-complementarity 2

# For systems with limited memory
./moducomp.py pipeline ./genomes ./output_lowmem --ncpus 8 --lowmem --calculate-complementarity 2
```

## Citation
If you use `moducomp` in your research, please cite:
Villada, JC. & Schulz, F. (2025). `moducomp`: A computational tool for metabolic module completeness of genomes and metabolic complementarity in microbiomes.
