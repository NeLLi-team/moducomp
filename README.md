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
- **Automatic resource monitoring** with timestamped logs tracking CPU usage, memory consumption, and runtime for reproducibility.
- **Consistent logging to stdout/stderr** with a per-command resource summary emitted at the end of each run.
- **Built-in validation (`moducomp validate`)** for scientific consistency checks across annotations, KO matrices, KPCT outputs, and complementarity reports.

## Installation (Recommended)

Make sure you have [`Pixi`](https://pixi.prefix.dev/latest/) installed:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

Install `moducomp` with `Pixi`:

```bash
pixi global install \
  -c conda-forge \
  -c bioconda \
  -c https://repo.prefix.dev/astrogenomics \
  moducomp
```

## Setup data (required)

`moducomp` needs the eggNOG-mapper database to run. The primary (recommended) way to download it is using the `download_eggnog_data.py` wrapper, which mirrors the upstream downloader behavior. For upstream details, see the eggNOG-mapper setup guide: [eggNOG-mapper database setup](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13#user-content-Setup).

```bash
export EGGNOG_DATA_DIR="/path/to/eggnog-data"
download_eggnog_data.py --eggnog-data-dir "$EGGNOG_DATA_DIR"
# equivalent:
# moducomp download-eggnog-data --eggnog-data-dir "$EGGNOG_DATA_DIR"
```

If `EGGNOG_DATA_DIR` is not set, the downloader defaults to `${XDG_DATA_HOME:-~/.local/share}/moducomp/eggnog`.

### Quick test

Small test data sets ship with `moducomp`. After installation you can confirm the pipeline by running:

```bash
moducomp test --ncpus 16 --calculate-complementarity 2 --eggnog-data-dir "$EGGNOG_DATA_DIR"
```

The test command runs in low-memory mode by default. If you have plenty of RAM and want full-memory mode, add `--fullmem` (or `--full-mem`).

### Developer install (Pixi)

If you want to download the code and develop locally:

```bash
git clone https://github.com/NeLLi-team/moducomp.git
cd moducomp
pixi install
```

## Quick check

```bash
moducomp --help
```

If you are running from the repository with `Pixi`:

```bash
pixi run python -m moducomp --help
```

You should see the command line help without errors.

## Usage

`moducomp` provides two main commands: `pipeline` and `analyze-ko-matrix`. You can run these commands using Pixi tasks defined in `pyproject.toml` or directly within the Pixi environment.

### Pipeline overview

The diagram below shows the main stages executed by ModuComp.

```mermaid
graph TD
    A([Start run]) --> B[Initialize logging and resource monitoring]
    B --> C{Input type}
    C -->|pipeline| D[Validate genome directory]
    C -->|analyze-ko-matrix| H[Load existing KO matrix]
    D --> E[Prepare genomes: adapt headers or copy to tmp]
    E --> F[Merge genomes into single FAA]
    F --> G[Run eggNOG-mapper if needed]
    G --> H[Create KO matrix kos_matrix.csv]
    H --> I[Convert KO matrix to KPCT input]
    I --> J[Run KPCT with parallel fallback]
    J --> K[Create module completeness matrix]
    K --> L{Complementarity requested?}
    L -->|Yes| M[Generate complementarity reports]
    L -->|No| N[Skip]
    M --> O[Write outputs and logs]
    N --> O
    O --> P[Optional cleanup of tmp]
    P --> Q([Pipeline complete])
```

### CLI options and defaults

This section lists all CLI options implemented today, along with their default values.

#### `pipeline` command (positional args: `genomedir`, `savedir`)

| Option | Default | Description |
| --- | --- | --- |
| `--ncpus`, `-n` | `16` | Number of CPU cores to use for eggNOG-mapper and KPCT. |
| `--calculate-complementarity`, `-c` | `0` | Complementarity size to compute (0 disables). |
| `--adapt-headers/--no-adapt-headers` | `false` | Adapt FASTA headers to `genome|protein_N`. |
| `--del-tmp/--keep-tmp` | `true` | Delete temporary files after completion. |
| `--lowmem/--fullmem` (`--low-mem/--full-mem`) | `fullmem` | Run eggNOG-mapper without `--dbmem` to reduce RAM. |
| `--verbose/--quiet` | `false` | Enable verbose progress output. |
| `--validate/--no-validate` | `validate` | Run post-run validation checks. |
| `--validate-report/--no-validate-report` | `validate-report` | Write `validation_report.json` in the output directory. |
| `--validate-strict/--validate-lenient` | `lenient` | Treat validation warnings as failures when strict. |
| `--log-level`, `-l` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`. |
| `--eggnog-data-dir` | `EGGNOG_DATA_DIR` | Path to eggNOG-mapper data (sets `EGGNOG_DATA_DIR`). |

#### `test` command (bundled test genomes)

| Option | Default | Description |
| --- | --- | --- |
| `--output-dir`, `-o` | `output_test_moducomp_<DATETIME>` | Output directory for test run. |
| `--ncpus`, `-n` | `2` | CPU cores for the test run. |
| `--calculate-complementarity`, `-c` | `2` | Complementarity size to compute (0 disables). |
| `--adapt-headers/--no-adapt-headers` | `false` | Adapt FASTA headers before the test. |
| `--del-tmp/--keep-tmp` | `true` | Delete temporary files after the test completes. |
| `--lowmem/--fullmem` (`--low-mem/--full-mem`) | `lowmem` | Low-memory mode is the default for tests. |
| `--verbose/--quiet` | `verbose` | Verbose output is the default for tests. |
| `--validate/--no-validate` | `validate` | Run post-run validation checks. |
| `--validate-report/--no-validate-report` | `validate-report` | Write `validation_report.json` in the output directory. |
| `--validate-strict/--validate-lenient` | `lenient` | Treat validation warnings as failures when strict. |
| `--log-level`, `-l` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`. |
| `--eggnog-data-dir` | `EGGNOG_DATA_DIR` | Path to eggNOG-mapper data (sets `EGGNOG_DATA_DIR`). |

#### `analyze-ko-matrix` command (positional args: `kos_matrix`, `savedir`)

| Option | Default | Description |
| --- | --- | --- |
| `--calculate-complementarity`, `-c` | `0` | Complementarity size to compute (0 disables). |
| `--kpct-outprefix` | `output_give_completeness` | Prefix for KPCT output files. |
| `--del-tmp/--keep-tmp` | `true` | Delete temporary files after completion. |
| `--ncpus`, `-n` | `16` | CPU cores for KPCT parallel processing. |
| `--verbose/--quiet` | `false` | Enable verbose progress output. |
| `--validate/--no-validate` | `validate` | Run post-run validation checks. |
| `--validate-report/--no-validate-report` | `validate-report` | Write `validation_report.json` in the output directory. |
| `--validate-strict/--validate-lenient` | `lenient` | Treat validation warnings as failures when strict. |
| `--log-level`, `-l` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

#### `validate` command (positional args: `savedir`)

| Option | Default | Description |
| --- | --- | --- |
| `--mode` | `auto` | Validation mode: `auto`, `pipeline`, or `ko-matrix`. |
| `--calculate-complementarity`, `-c` | `auto-detect` | Expected complementarity size (0 disables). |
| `--kpct-outprefix` | `output_give_completeness` | KPCT output prefix used during analysis. |
| `--strict/--lenient` | `lenient` | Treat warnings as failures when strict. |
| `--report` | _none_ | Write JSON validation report to this path. |
| `--verbose/--quiet` | `false` | Enable verbose progress output. |
| `--log-level`, `-l` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`. |

#### `download-eggnog-data` command

| Option | Default | Description |
| --- | --- | --- |
| `--eggnog-data-dir` | `${XDG_DATA_HOME:-~/.local/share}/moducomp/eggnog` | Destination for eggNOG-mapper data (sets `EGGNOG_DATA_DIR`). |
| `--log-level`, `-l` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR`. |
| `--verbose/--quiet` | `verbose` | Stream downloader output to the console. |

### Performance and parallel processing

`moducomp` includes **parallel processing capabilities** for the KPCT (KEGG Pathways Completeness Tool) analysis, which can significantly improve performance for large datasets:

- **Automatic chunking**: Input files are automatically split into chunks for parallel processing
- **Checkpointing**: Resume capability if processing is interrupted - already processed chunks are automatically detected and skipped
- **Fallback mechanism**: If parallel processing fails, the system automatically falls back to sequential processing
- **Configurable CPU usage**: Use the `--ncpus` parameter to control how many CPU cores to use

**CPU Configuration**:
- The `--ncpus` parameter controls the number of CPU cores used for both eggNOG-mapper annotation and KPCT analysis
- For KPCT parallel processing, the system creates the same number of chunks as CPU cores specified
- Example: `--ncpus 8` will use 8 cores and create 8 chunks for optimal parallel processing

### Validation (QC)

Use the built-in validator to check scientific consistency across outputs after a run. The validator compares:
- KO sets and counts between eggNOG-mapper annotations and `kos_matrix.csv`
- KO sets between `kos_matrix.csv` and `ko_file_for_kpct.txt`
- KPCT contigs vs pathways outputs
- Module completeness ranges and combination naming
- Complementarity reports versus module completeness values
- Protein provenance fields (pipeline mode) or placeholders (KO-matrix mode)

Example:

```bash
# Validation runs by default after pipeline/analyze/test.
# Use --no-validate to disable or --no-validate-report to skip JSON output.
# When validation reports errors (or warnings in strict mode), the command exits non-zero.

# Validate a pipeline run and write a JSON report
moducomp validate /path/to/output --mode pipeline --report /path/to/output/validation_report.json

# Validate KO-matrix mode outputs (non-default KPCT prefix)
moducomp validate /path/to/output --mode ko-matrix --kpct-outprefix my_prefix

# Treat warnings as failures
moducomp validate /path/to/output --strict
```

### ⚠️ Important note 1

**Prepare FAA files**: Ensure FAA headers are in the form `>genomeName|proteinId`, or use the `--adapt-headers` option to format your headers into `>fileName_prefix|protein_id_counter`.

### ⚠️ Important note 2

`moducomp` is specifically designed for large scale analysis of microbiomes with hundreds of members, and works on Linux systems with at least **64GB of RAM**. Nevertheless, it can be run on **smaller systems with less RAM, using the flag `--lowmem` (`--low-mem`) when running the `pipeline` command**. The `test` command uses low-memory mode by default and can be switched to full memory with `--fullmem` (`--full-mem`).

### Notes on bundled test data

You can override the bundled data location with `MODUCOMP_DATA_DIR`.
When working from source, the bundled test genomes live at `moducomp/data/test_genomes`.

`download_eggnog_data.py` is exposed by `moducomp` as a convenience wrapper for the eggnog-mapper downloader and is available in the Pixi environment (including `pixi global` installs).

Pixi task (supports passing a custom location):

```bash
export EGGNOG_DATA_DIR=/path/to/eggnog-data
pixi run download-eggnog-data --eggnog-data-dir /path/to/eggnog-data
```

Global install shortcut (also supports `--eggnog-data-dir`):

```bash
download-eggnog-data --eggnog-data-dir /path/to/eggnog-data
```

### Running with your samples

If you are running from the repository with `Pixi`, replace `moducomp` below with `pixi run python -m moducomp`.

#### `pipeline` command

Use the `pipeline` command to process a directory of genome FAA files from scratch.

```bash
moducomp pipeline \
    /path/to/your/faa_files \
    /path/to/your/output_directory \
    --ncpus <number_of_cpus_to_use> \
    --calculate-complementarity <N>  # 0 to disable, 2 for 2-member, 3 for 3-member complementarity.
    # Optional flags:
    # --lowmem/--fullmem          # Optional: Use low-mem if you have less than 64GB of RAM (default is full mem)
    # --adapt-headers             # If your FASTA headers need modification
    # --del-tmp/--keep-tmp        # Delete or keep temporary files
    # --eggnog-data-dir /path     # If EGGNOG_DATA_DIR is not set
    # --verbose                   # Enable verbose output with detailed progress information
```

#### `analyze-ko-matrix` command

Use the `analyze-ko-matrix` command if you already have a KO matrix file (CSV format, where rows are genomes/combinations, columns are KOs, and values are KO counts).

The KO matrix file should have a `taxon_oid` column for genome identifiers, and subsequent columns for each KO (e.g., `K00001`, `K00002`) with integer counts.

```bash
moducomp analyze-ko-matrix \
    /path/to/your/kos_matrix.csv \
    /path/to/your/output_directory \
    --ncpus <number_of_cpus_to_use> \
    --calculate-complementarity <N>  # 0 to disable, 2 for 2-member, 3 for 3-member complementarity.

    # Optional flags:
    # --keep-tmp                  # Keep temporary files
    # --verbose                   # Enable verbose output with detailed progress information
```

### Parallel processing features

`moducomp` includes advanced parallel processing capabilities for improved performance:

#### KPCT parallel processing

When using the `--ncpus` parameter with a value greater than 1, `moducomp` automatically enables parallel processing for the KPCT (KEGG Pathways Completeness Tool) analysis:

- **Automatic chunking**: Input files are split into `ncpus` chunks for optimal load balancing
- **Concurrent processing**: Multiple chunks are processed simultaneously using `multiprocessing.ProcessPoolExecutor`
- **Resume capability**: If processing is interrupted, completed chunks are automatically detected and skipped on restart
- **Automatic fallback**: If parallel processing fails, the system seamlessly falls back to sequential processing

#### Performance tips

- **CPU cores**: Start with `--ncpus 8` for moderate datasets, increase to `--ncpus 16` or higher for large datasets
- **Memory considerations**: Each parallel worker requires memory; reduce `--ncpus` if you encounter memory issues
- **Large datasets**: For datasets with hundreds of genomes, parallel processing can reduce KPCT analysis time by 50-80%

#### Example with parallel processing

```bash
# For large datasets with sufficient resources
moducomp pipeline ./large_genome_collection ./output_large --ncpus 32 --calculate-complementarity 3

# For moderate datasets with verbose output
moducomp analyze-ko-matrix ./ko_matrix.csv ./output_moderate --ncpus 16 --calculate-complementarity 2 --verbose

# For systems with limited memory
moducomp pipeline ./genomes ./output_lowmem --ncpus 8 --lowmem --calculate-complementarity 2
```

## Expected outputs

The sections below describe the expected output files, naming conventions, and the column-level meaning of each file. These details are the same for `moducomp pipeline` and `moducomp test` (pipeline mode), and the subset noted for `moducomp analyze-ko-matrix` (KO-matrix mode).

**Naming conventions**

Genome identifiers are stored as `taxon_oid`. In pipeline mode, ModuComp expects protein headers in the format `genome_id|protein_id`. If you set `--adapt-headers`, ModuComp rewrites headers to `>genomeName|protein_N`, where `genomeName` is the FAA filename stem. Combination identifiers use `__` (double underscore), for example `GenomeA__GenomeB`, and `n_members` in `module_completeness.tsv` records the size of each combination.

**Pipeline mode outputs (`moducomp pipeline`, `moducomp test`)**

- `emapper_out.emapper.annotations`: Full eggNOG-mapper annotations. The `#query` column must match `genome_id|protein_id`. `KEGG_ko` entries are prefixed `ko:KXXXXX` and are converted to `KXXXXX` for downstream matrices.
- `kos_matrix.csv`: Genome × KO count matrix. Columns: `taxon_oid` followed by KO IDs (e.g., `K00001`). Values are integer protein counts per KO.
- `ko_file_for_kpct.txt`: KPCT input file. Each line starts with `taxon_oid` followed by the set of KO IDs present in that genome or combination. If `--calculate-complementarity` is `N>=2`, combinations up to `N` are included as `GenomeA__GenomeB`.
- `output_give_completeness_contigs.with_weights.tsv`: KPCT module results per genome/combination. Columns: `contig` (genome/combination ID), `module_accession`, `completeness` (0–100), `pathway_name`, `pathway_class`, `matching_ko` (KO weights), `missing_ko`.
- `output_give_completeness_pathways.with_weights.tsv`: Same rows and order as the contigs file, but without the `contig` column. This is provided for compatibility with legacy tools; prefer the contigs file when you need genome-level provenance.
- `module_completeness.tsv`: Pivoted module completeness matrix. Columns: `n_members`, `taxon_oid`, followed by KEGG module IDs (`M00001`, …). Values are numeric percentages in the range 0–100.
- `module_completeness_complementarity_Nmember.tsv`: Complementarity report for `N`-member combinations (only when `--calculate-complementarity N` is set). Columns: `taxon_oid_1..N`, `completeness_taxon_oid_1..N`, `module_id`, `module_name`, `pathway_class`, `matching_ko`, `proteins_taxon_oid_1..N`. Protein fields list contributing proteins per KO (from eggNOG-mapper) as `{'KXXXXX': 'genome|protein'}`.
- `logs/moducomp.log`: Detailed run log with structured progress messages and per-command resource summaries.
- `logs/resource_usage_YYYYMMDD_HHMMSS.log`: Resource monitoring log capturing wall time, CPU time, CPU utilization, peak RAM, and exit code for each monitored command.
- `tmp/` (only if `--keep-tmp`): Intermediate files such as `merged_genomes.faa`, `emapper_output/`, and KPCT chunk outputs.
- `validation_report.json` (default when validation is enabled): JSON report produced by the validator.

**KO-matrix mode outputs (`moducomp analyze-ko-matrix`)**

- `kos_matrix.csv`: A copy of the input KO matrix (same format as above).
- `ko_file_for_kpct.txt`: KPCT input generated from the KO matrix. If `--calculate-complementarity` is set, combination lines are added using `GenomeA__GenomeB` identifiers.
- `output_give_completeness_contigs.with_weights.tsv`: KPCT module results per genome/combination (same format as pipeline mode).
- `output_give_completeness_pathways.with_weights.tsv`: Same rows as the contigs file, without the `contig` column.
- `module_completeness.tsv`: Module completeness matrix (same format as pipeline mode).
- `module_completeness_complementarity_Nmember.tsv`: Complementarity report. Protein contribution columns are filled with `No protein data available for <genome>` because no eggNOG-mapper annotations are available in KO-matrix mode.
- `logs/moducomp.log` and `logs/resource_usage_YYYYMMDD_HHMMSS.log`: Standard run logs and resource summaries.
- `validation_report.json` (default when validation is enabled): JSON report produced by the validator.

## Citation
Villada, JC. & Schulz, F. (2025). Assessment of metabolic module completeness of genomes and metabolic complementarity in microbiomes with `moducomp` . `moducomp` (v0.5.1) Zenodo. https://doi.org/10.5281/zenodo.16116092
