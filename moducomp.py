#!/usr/bin/env python
"""
moducomp: Metabolic module completeness of genomes and metabolic complementarity in microbiomes

This module provides a comprehensive bioinformatics pipeline for analyzing metabolic
module completeness in microbial genomes and identifying complementarity patterns
in microbial communities.

Key Features:
- Protein annotation using eggNOG-mapper to obtain KO (KEGG Orthology) terms
- Mapping of KO terms to KEGG metabolic modules using KPCT
- Parallel processing support for improved performance
- Module completeness analysis for individual genomes
- Complementarity analysis for N-member genome combinations
- Protein-level tracking for module completion

The tool supports two main workflows:
1. Complete pipeline: FAA files → annotations → KO matrix → module analysis
2. KO matrix analysis: Pre-existing KO matrix → module analysis

Author: Juan C. Villada - US DOE Joint Genome Institute - Lawrence Berkeley National Lab
License: See LICENSE.txt
Version: See pixi.toml for current version
"""

import datetime
import glob
import itertools
import logging
import multiprocessing
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from sys import exit
from typing import Any, Dict, List, Optional, Set, Tuple, Union


import pandas as pd
import typer

def conditional_output(message: str, color: str = "white", verbose: bool = True) -> None:
    """
    Print message to terminal only if verbose mode is enabled.

    Parameters
    ----------
    message : str
        Message to display
    color : str, optional
        Color for the message, by default "white"
    verbose : bool, optional
        Whether to display the message, by default True
    """
    if verbose:
        typer.secho(message, fg=color)

def run_subprocess_with_logging(
    cmd: List[str],
    logger: Optional[logging.Logger] = None,
    description: str = "Running command",
    verbose: bool = True
) -> Tuple[int, str, str]:
    """
    Run a subprocess command with proper stdout/stderr streaming to both console and log file.

    Parameters
    ----------
    cmd : List[str]
        Command and arguments to run
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress
    description : str, optional
        Description of what the command does for user display

    Returns
    -------
    Tuple[int, str, str]
        Return code, stdout, and stderr as strings
    """
    if logger:
        logger.info(f"{description}: {' '.join(cmd)}")

    # Only show detailed command info in verbose mode
    conditional_output(f"🔧 {description}", "yellow", verbose)
    conditional_output(f"   Command: {' '.join(cmd)}", "blue", verbose)

    if verbose:
        print(f"[DEBUG] Starting subprocess with command: {cmd[0]}")
        print(f"[DEBUG] Working directory: {os.getcwd()}")
        sys.stdout.flush()

    try:
        def stream_reader(stream, q, stream_type):
            """Read from stream and put lines in queue"""
            try:
                while True:
                    line = stream.readline()
                    if not line:
                        break
                    line = line.rstrip('\n\r')
                    q.put((stream_type, line))  # Put all lines, even empty ones
                stream.close()
            except Exception:
                pass

        # Start the process
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=0,  # Unbuffered
            universal_newlines=True
        )

        # Create queues for stdout and stderr
        stdout_queue = queue.Queue()
        stderr_queue = queue.Queue()

        # Start threads to read stdout and stderr
        stdout_thread = threading.Thread(
            target=stream_reader,
            args=(process.stdout, stdout_queue, 'stdout')
        )
        stderr_thread = threading.Thread(
            target=stream_reader,
            args=(process.stderr, stderr_queue, 'stderr')
        )

        stdout_thread.daemon = True
        stderr_thread.daemon = True
        stdout_thread.start()
        stderr_thread.start()

        stdout_lines = []
        stderr_lines = []
        last_output_time = time.time()
        progress_interval = 1800  # Show progress every 30 minutes

        if verbose:
            print(f"[DEBUG] Starting to monitor subprocess output...")
            print(f"[DEBUG] Process PID: {process.pid}")
            print(f"[DEBUG] Process poll status: {process.poll()}")
            sys.stdout.flush()

        # Read from queues and process output in real-time
        while process.poll() is None or not stdout_queue.empty() or not stderr_queue.empty():
            current_time = time.time()
            output_received = False

            # Check stdout queue
            try:
                stream_type, line = stdout_queue.get_nowait()
                if stream_type == 'stdout':
                    stdout_lines.append(line)
                    # Stream to console immediately
                    if verbose:
                        print(line, flush=True)
                        if len(stdout_lines) <= 3:  # Debug first few lines
                            print(f"[DEBUG] Received stdout line {len(stdout_lines)}: {repr(line)}", flush=True)
                    if logger:
                        logger.info(line)
                    last_output_time = current_time
                    output_received = True
            except queue.Empty:
                pass

            # Check stderr queue
            try:
                stream_type, line = stderr_queue.get_nowait()
                if stream_type == 'stderr':
                    stderr_lines.append(line)
                    # Stream to console immediately
                    if verbose:
                        print(f"[STDERR] {line}", file=sys.stderr, flush=True)
                        if len(stderr_lines) <= 3:  # Debug first few stderr lines
                            print(f"[DEBUG] Received stderr line {len(stderr_lines)}: {repr(line)}", flush=True)
                    if logger:
                        logger.warning(line)
                    last_output_time = current_time
                    output_received = True
            except queue.Empty:
                pass

            # Show progress message if no output for a while
            if not output_received and current_time - last_output_time > progress_interval:
                if verbose:
                    elapsed = int(current_time - last_output_time)
                    print(f"   ... still running (no output for {elapsed}s)", flush=True)
                    print(f"[DEBUG] Process status: {process.poll()}, PID: {process.pid}", flush=True)
                if logger:
                    logger.info(f"Process still running, no output for {int(current_time - last_output_time)} seconds")
                last_output_time = current_time

            # Small delay to prevent busy waiting
            time.sleep(0.05)  # Reduced delay for more responsive output

        # Wait for threads to complete
        stdout_thread.join(timeout=1.0)
        stderr_thread.join(timeout=1.0)

        # Get any remaining items from queues
        while not stdout_queue.empty():
            try:
                stream_type, line = stdout_queue.get_nowait()
                if stream_type == 'stdout':
                    stdout_lines.append(line)
                    if verbose:
                        print(line, flush=True)
                    if logger:
                        logger.info(line)
            except queue.Empty:
                break

        while not stderr_queue.empty():
            try:
                stream_type, line = stderr_queue.get_nowait()
                if stream_type == 'stderr':
                    stderr_lines.append(line)
                    if verbose:
                        print(line, file=sys.stderr, flush=True)
                    if logger:
                        logger.warning(line)
            except queue.Empty:
                break

        returncode = process.returncode
        stdout_str = '\n'.join(stdout_lines)
        stderr_str = '\n'.join(stderr_lines)

        if logger:
            logger.info(f"Command completed with return code: {returncode}")

        return returncode, stdout_str, stderr_str

    except Exception as e:
        error_msg = f"Exception running command {' '.join(cmd)}: {str(e)}"
        if logger:
            logger.error(error_msg)
        print(f"[ERROR] {error_msg}", file=sys.stderr)
        return -1, "", str(e)


def setup_resource_logging(savedir: str) -> str:
    """
    Set up resource usage logging file with timestamp.

    Parameters
    ----------
    savedir : str
        Directory to save the resource log file

    Returns
    -------
    str
        Path to the resource log file
    """
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    resource_log_file = os.path.join(savedir, f"resource_usage_{timestamp}.log")

    # Create header for resource log file
    with open(resource_log_file, 'w') as f:
        f.write("# ModuComp Resource Usage Report\n")
        f.write(f"# Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# Format: Timestamp | Command | Elapsed_Time(s) | User_Time(s) | System_Time(s) | CPU_Percent | Max_RAM(GB) | Exit_Code\n")
        f.write("#" + "="*120 + "\n\n")

    return resource_log_file


def run_subprocess_with_resource_monitoring(
    cmd: List[str],
    resource_log_file: str,
    logger: Optional[logging.Logger] = None,
    description: str = "Running command",
    verbose: bool = True
) -> Tuple[int, str, str]:
    """
    Run a subprocess command with resource monitoring using /usr/bin/time.

    Parameters
    ----------
    cmd : List[str]
        Command and arguments to run
    resource_log_file : str
        Path to the resource usage log file
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress
    description : str, optional
        Description of what the command does for user display
    verbose : bool, optional
        Whether to display detailed output

    Returns
    -------
    Tuple[int, str, str]
        Return code, stdout, and stderr as strings
    """
    # Create time format string for detailed resource monitoring
    time_format = "Command: %C\\nElapsed Time (wall clock): %E\\nElapsed Time (seconds): %e\\nUser Time (seconds): %U\\nSystem Time (seconds): %S\\nCPU Percentage: %P\\nMaximum RAM (KB): %M\\nExit Status: %x\\nMajor Page Faults: %F\\nMinor Page Faults: %R\\nContext Switches (voluntary): %w\\nContext Switches (involuntary): %c\\n"

    # Create temporary file for time output
    temp_time_file = f"{resource_log_file}.tmp"

    # Wrap the original command with /usr/bin/time
    time_cmd = [
        "/usr/bin/time",
        "-f", time_format,
        "-o", temp_time_file
    ] + cmd

    if logger:
        logger.info(f"{description}: {' '.join(cmd)}")
        logger.info(f"Resource monitoring enabled, output will be logged to: {resource_log_file}")

    conditional_output(f"🔧 {description} (with resource monitoring)", "yellow", verbose)
    conditional_output(f"   Command: {' '.join(cmd)}", "blue", verbose)

    start_time = datetime.datetime.now()

    # Run the command with the existing subprocess function but with time wrapper
    # Note: Pass None as logger to avoid logging the verbose time command
    returncode, stdout, stderr = run_subprocess_with_logging(
        time_cmd, None, description, verbose
    )

    end_time = datetime.datetime.now()

    # Parse resource usage from temporary file
    resource_info = {}
    if os.path.exists(temp_time_file):
        try:
            with open(temp_time_file, 'r') as f:
                time_output = f.read()

            # Parse the time output
            for line in time_output.strip().split('\n'):
                if ':' in line:
                    key, value = line.split(':', 1)
                    resource_info[key.strip()] = value.strip()

            # Extract key metrics for summary
            elapsed_seconds = resource_info.get('Elapsed Time (seconds)', 'N/A')
            user_time = resource_info.get('User Time (seconds)', 'N/A')
            system_time = resource_info.get('System Time (seconds)', 'N/A')
            cpu_percent = resource_info.get('CPU Percentage', 'N/A')
            max_ram_kb = resource_info.get('Maximum RAM (KB)', 'N/A')
            exit_status = resource_info.get('Exit Status', str(returncode))

            # Convert RAM from KB to GB
            try:
                max_ram_gb = float(max_ram_kb) / (1024 * 1024) if max_ram_kb != 'N/A' else 'N/A'
                max_ram_gb_str = f"{max_ram_gb:.3f}" if max_ram_gb != 'N/A' else 'N/A'
            except (ValueError, TypeError):
                max_ram_gb_str = 'N/A'

            # Write summary to resource log
            timestamp_str = start_time.strftime("%Y-%m-%d %H:%M:%S")
            cmd_str = ' '.join(cmd)

            with open(resource_log_file, 'a') as f:
                f.write(f"{timestamp_str} | {cmd_str} | {elapsed_seconds} | {user_time} | {system_time} | {cpu_percent} | {max_ram_gb_str} | {exit_status}\n")
                f.write("# Detailed resource information:\n")
                for key, value in resource_info.items():
                    f.write(f"# {key}: {value}\n")
                f.write("\n")

            # Display resource summary
            if verbose:
                conditional_output("\n📊 Resource Usage Summary:", "cyan", verbose)
                conditional_output(f"   ⏱️  Wall Clock Time: {elapsed_seconds}s", "white", verbose)
                conditional_output(f"   💻 CPU Time (User): {user_time}s", "white", verbose)
                conditional_output(f"   🔧 CPU Time (System): {system_time}s", "white", verbose)
                conditional_output(f"   📈 CPU Usage: {cpu_percent}", "white", verbose)
                conditional_output(f"   🧠 Peak RAM Usage: {max_ram_gb_str} GB", "white", verbose)
                conditional_output(f"   ✅ Exit Code: {exit_status}\n", "white", verbose)

            if logger:
                logger.info(f"Resource usage - Wall time: {elapsed_seconds}s, CPU: {cpu_percent}, Peak RAM: {max_ram_gb_str} GB")

        except Exception as e:
            if logger:
                logger.warning(f"Failed to parse resource usage: {str(e)}")
            conditional_output(f"⚠️  Warning: Could not parse resource usage: {str(e)}", "yellow", verbose)

        # Clean up temporary file
        try:
            os.remove(temp_time_file)
        except OSError:
            pass

    else:
        if logger:
            logger.warning("Resource monitoring file not found")
        conditional_output("⚠️  Warning: Resource monitoring output not found", "yellow", verbose)

    return returncode, stdout, stderr


def log_final_resource_summary(resource_log_file: str, total_start_time: float, logger: Optional[logging.Logger] = None, verbose: bool = True) -> None:
    """
    Log final resource summary at the end of the pipeline.

    Parameters
    ----------
    resource_log_file : str
        Path to the resource usage log file
    total_start_time : float
        Start time of the entire pipeline (from time.time())
    logger : Optional[logging.Logger], optional
        Logger instance
    verbose : bool, optional
        Whether to display summary to console
    """
    total_elapsed = time.time() - total_start_time
    end_time = datetime.datetime.now()

    with open(resource_log_file, 'a') as f:
        f.write("\n" + "="*120 + "\n")
        f.write("# PIPELINE SUMMARY\n")
        f.write(f"# Pipeline completed at: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Total pipeline elapsed time: {total_elapsed:.2f} seconds ({total_elapsed/60:.2f} minutes)\n")
        f.write("="*120 + "\n")

    if verbose:
        conditional_output("\n🎯 Resource Usage Summary Logged", "green", verbose)
        conditional_output(f"   📋 Full report saved to: {resource_log_file}", "cyan", verbose)
        conditional_output(f"   ⏱️  Total pipeline time: {total_elapsed:.2f}s ({total_elapsed/60:.2f}min)", "white", verbose)

    if logger:
        logger.info(f"Resource usage summary completed. Total time: {total_elapsed:.2f}s")
        logger.info(f"Resource log saved to: {resource_log_file}")


def parse_emapper_annotations(emapper_file_path: str, logger: Optional[logging.Logger] = None) -> Dict[str, Dict[str, List[str]]]:
    """
    Parses an emapper annotation file to extract protein to KO mappings.

    Parameters
    ----------
    emapper_file_path : str
        Path to the emapper.annotations file.
    logger : Optional[logging.Logger]
        Logger instance.

    Returns
    -------
    Dict[str, Dict[str, List[str]]]
        A nested dictionary mapping:
        genome_id -> {ko_id -> [protein_id1, protein_id2, ...]}
        This allows tracking which proteins in each genome map to which KO terms.
    """
    genome_ko_proteins: Dict[str, Dict[str, List[str]]] = {}

    if logger:
        logger.info(f"Parsing emapper annotations from: {emapper_file_path}")

    try:

        emapper_df = pd.read_csv(
            emapper_file_path,
            sep="\t",
            skiprows=4,
            skipfooter=3,
            engine='python',
            usecols=["#query", "KEGG_ko"]
        )


        emapper_df = emapper_df[emapper_df["KEGG_ko"] != "-"]
        emapper_df = emapper_df.dropna(subset=["KEGG_ko", "#query"])

        if logger:
            logger.info(f"Processing {len(emapper_df)} rows with KO annotations")

        for _, row in emapper_df.iterrows():
            protein_id = str(row["#query"])
            kegg_ko_str = str(row["KEGG_ko"])


            parts = protein_id.split('|', 1)
            if len(parts) != 2:

                if logger:
                    logger.warning(f"Protein ID doesn't follow expected format (genome_id|protein_id): {protein_id}")
                continue

            genome_id = parts[0]


            if genome_id not in genome_ko_proteins:
                genome_ko_proteins[genome_id] = {}


            if kegg_ko_str and kegg_ko_str != '-':
                kos = kegg_ko_str.split(',')
                for ko in kos:

                    ko_cleaned = ko.replace('ko:', '').strip()


                    if '(' in ko_cleaned:
                        ko_cleaned = ko_cleaned.split('(')[0].strip()

                    if ko_cleaned:

                        if ko_cleaned not in genome_ko_proteins[genome_id]:
                            genome_ko_proteins[genome_id][ko_cleaned] = []

                        genome_ko_proteins[genome_id][ko_cleaned].append(protein_id)

        if logger:
            total_genomes = len(genome_ko_proteins)
            total_kos = sum(len(kos) for kos in genome_ko_proteins.values())
            logger.info(f"Successfully parsed protein mappings for {total_genomes} genomes with {total_kos} total KO annotations")


            if total_genomes > 0:
                sample_genome = next(iter(genome_ko_proteins))
                sample_kos = list(genome_ko_proteins[sample_genome].keys())[:5]
                logger.debug(f"Sample genome {sample_genome} has KOs: {sample_kos}")


                key_kos = ['K00878', 'K00941']
                for key_ko in key_kos:
                    for genome, ko_dict in genome_ko_proteins.items():
                        if key_ko in ko_dict:
                            logger.debug(f"Found {key_ko} in genome {genome} with {len(ko_dict[key_ko])} proteins")
                            break

    except FileNotFoundError:
        if logger:
            logger.error(f"Emapper annotation file not found: {emapper_file_path}")
        return {}
    except pd.errors.EmptyDataError:
        if logger:
            logger.error(f"Emapper annotation file is empty or resulted in an empty DataFrame after skipping rows/footer: {emapper_file_path}")
        return {}
    except ValueError as ve:
        if logger:
            logger.error(f"ValueError during pandas parsing for {emapper_file_path}: {ve}")
        return {}
    except Exception as e:
        if logger:
            logger.error(f"An unexpected error occurred while parsing {emapper_file_path}: {e}")
        return {}

    return genome_ko_proteins

def get_module_protein_contributions(
    module_id: str,
    module_kos: List[str],
    genome_ko_proteins: Dict[str, Dict[str, List[str]]],
    logger: Optional[logging.Logger] = None
) -> Dict[str, Dict[str, List[str]]]:
    """
    Get the protein contributions for a module from each genome.

    Parameters
    ----------
    module_id : str
        The ID of the module (e.g., "M00001").
    module_kos : List[str]
        List of KO IDs that belong to this module.
    genome_ko_proteins : Dict[str, Dict[str, List[str]]]
        The nested dictionary mapping genome_id -> {ko_id -> [protein_ids]}.
    logger : Optional[logging.Logger], optional
        Logger instance.

    Returns
    -------
    Dict[str, Dict[str, List[str]]]
        A nested dictionary mapping:
        genome_id -> {ko_id -> [protein_ids]}
        Only including KOs that are part of the module.
    """
    result = {}

    for genome_id, ko_proteins in genome_ko_proteins.items():
        genome_contributions = {}

        for ko in module_kos:
            if ko in ko_proteins:
                genome_contributions[ko] = ko_proteins[ko]

        if genome_contributions:
            result[genome_id] = genome_contributions

    if logger:
        contributing_genomes = len(result)
        if contributing_genomes > 0:
            total_proteins = sum(sum(len(proteins) for proteins in kos.values()) for kos in result.values())
            logger.debug(f"Found {contributing_genomes} genomes contributing {total_proteins} proteins to module {module_id}")

    return result

def format_protein_contributions(genome_id: str, contributions: Dict[str, List[str]]) -> str:
    """
    Format protein contributions for a single genome into a string.

    Parameters
    ----------
    genome_id : str
        The ID of the genome.
    contributions : Dict[str, List[str]]
        Dictionary mapping KO IDs to lists of protein IDs.

    Returns
    -------
    str
        Formatted string representing the contributions.
        Format: {'KO1': 'protein1, protein2', 'KO2': 'protein3, protein4'}
    """
    formatted_contributions = {}

    for ko, proteins in contributions.items():

        formatted_contributions[ko] = ', '.join(proteins)


    if formatted_contributions:
        return str(formatted_contributions)
    else:
        return f"No proteins from {genome_id} contribute to this module"



def setup_logging(savedir: str) -> logging.Logger:
    """
    Set up logging to both console and file.

    Parameters
    ----------
    savedir : str
        Directory where log file will be saved

    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    os.makedirs(savedir, exist_ok=True)

    log_file = os.path.join(savedir, f"moducomp_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    log_file_dir = os.path.dirname(log_file)

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir, exist_ok=True)


    logger = logging.getLogger('moducomp')
    logger.setLevel(logging.DEBUG)


    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')


    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)


    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)


    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    logger.info(f"Log file created at: {log_file}")
    return logger

def greetings(verbose: bool = True):
    """
    Print welcome message for the ModuComp application.

    This function displays a formatted greeting message to the user
    when the application starts.
    """
    my_name = "\n📊 ModuComp"
    conditional_output(my_name, "green", verbose)


def get_path_to_each_genome(genomedir: str):
    """
    Get paths to all FAA files in the genome directory.

    Parameters
    ----------
    genomedir : str
        Directory containing genome files in FAA format

    Returns
    -------
    List[str]
        List of absolute paths to FAA files
    """
    return glob.glob(pathname=f"{genomedir}/*.faa")


def how_many_genomes(genomedir: str, verbose: bool = True):
    """
    Count and display the number of FAA files in the genome directory.

    This function checks for .faa files in the specified directory and
    displays a success or error message to the user. If no files are found,
    the program exits.

    Parameters
    ----------
    genomedir : str
        Directory containing genome files in FAA format
    verbose : bool
        Whether to display detailed output
    """
    n_files = len(get_path_to_each_genome(genomedir))
    if n_files > 0:
        conditional_output(f"✅ [OK] {n_files} faa files were found in '{genomedir}'\n", "white", verbose)
    else:
        # Always show errors regardless of verbose setting
        typer.secho(f"❌ [ERROR] No FAA files were found in '{genomedir}'\n", fg="red")
        exit()


def create_output_dir(savedir: str, verbose: bool = True):
    """
    Create the output directory if it doesn't exist.

    Parameters
    ----------
    savedir : str
        Path to the output directory to create
    verbose : bool
        Whether to display detailed output
    """
    conditional_output("\n📊 Creating output directory", "green", verbose)
    if os.path.exists(savedir):
        conditional_output(f"✅ [OK] Output directory already exists at: {savedir}\n", "white", verbose)
    else:
        os.makedirs(savedir, exist_ok=True)
        conditional_output(f"✅ [OK] Output directory created at: {savedir}\n", "white", verbose)


def get_tmp_dir(savedir:str) -> str:
    """
    Get the path to the temporary directory.

    Parameters
    ----------
    savedir : str
        Base output directory

    Returns
    -------
    str
        Path to the temporary directory within savedir
    """
    tmp_dir_path = f"{savedir}/tmp"
    return(tmp_dir_path)


def create_tmp_dir(savedir: str, verbose: bool = True):
    """
    Create the temporary directory if it doesn't exist.

    Parameters
    ----------
    savedir : str
        Base output directory where temporary directory will be created
    verbose : bool
        Whether to display detailed output
    """
    conditional_output("\n📊 Creating tmp dir", "green", verbose)
    tmp_dir_path = get_tmp_dir(savedir)
    if (os.path.exists(tmp_dir_path)):
        conditional_output(f"✅ [OK] Tmp directory already exists at: {tmp_dir_path}\n", "white", verbose)
    else:
        os.mkdir(tmp_dir_path)
        conditional_output(f"✅ [OK] Tmp directory created at: {tmp_dir_path}\n", "white", verbose)


def adapt_fasta_headers(genomedir: str, savedir: str, verbose: bool = True) -> None:
    """
    Modify FASTA headers to follow a consistent naming convention.

    This function processes all FAA files in the genome directory and modifies
    their headers to follow the format: >genomeName|protein_N where N is a
    sequential counter. This is required for proper downstream processing
    by eggNOG-mapper.

    Parameters
    ----------
    genomedir : str
        Directory containing original genome files in FAA format
    savedir : str
        Output directory where modified files will be saved
    verbose : bool
        Whether to display detailed output
    """
    conditional_output("\n📊 Modifying fasta headers", "green", verbose)
    path_to_each_genome = get_path_to_each_genome(genomedir)
    output_dir = f"{get_tmp_dir(savedir)}/faa"
    if os.path.exists(output_dir):
        conditional_output(f"✅ [OK] Fasta headers already modified at: {output_dir}\n", "white", verbose)
        return

    os.mkdir(output_dir)
    conditional_output("Processing fasta files and modifying headers...", "yellow", verbose)
    for each_file in path_to_each_genome:
        with open(each_file) as infile:
            with open(f"{output_dir}/{os.path.basename(each_file)}", "w") as outfile:
                genome_name_for_header = os.path.basename(each_file).split(".")[0]
                i=1
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(f">{genome_name_for_header}|protein_{i}\n")
                        i+=1
                    else:
                        outfile.write(line)
    conditional_output(f"✅ [OK] Fasta headers modified at: {output_dir}\n", "white", verbose)


def copy_faa_to_tmp(genomedir: str, savedir: str, verbose: bool = True) -> None:
    """
    Copy FAA files from the genome directory to the temporary directory.

    This function copies all .faa files from the input genome directory
    to the temporary processing directory without modification. Used when
    header adaptation is not required.

    Parameters
    ----------
    genomedir : str
        Source directory containing genome files in FAA format
    savedir : str
        Output directory where files will be copied to tmp/faa subdirectory
    verbose : bool
        Whether to display detailed output
    """
    conditional_output("\n📊 Copying faa files to tmp dir", "green", verbose)
    path_to_each_genome = get_path_to_each_genome(genomedir)
    output_dir = f"{get_tmp_dir(savedir)}/faa"
    if os.path.exists(output_dir):
        conditional_output(f"✅ [OK] Fasta files already exist at: {output_dir}\n", "white", verbose)
        return

    os.mkdir(output_dir)
    conditional_output("Copying genome files to temporary directory...", "yellow", verbose)
    for each_file in path_to_each_genome:
        shutil.copy(each_file, output_dir)
    conditional_output(f"✅ [OK] Fasta files copied to: {output_dir}\n", "white", verbose)


def merge_genomes(savedir: str, logger: Optional[logging.Logger] = None, verbose: bool = True) -> bool:
    """
    Merge individual genome FAA files into a single file.

    Parameters
    ----------
    savedir : str
        Directory where outputs will be saved
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress
    verbose : bool
        Whether to display detailed output

    Returns
    -------
    bool
        True if the merged file was created or already exists, False otherwise
    """
    conditional_output("\n📊 Merging genomes", "green", verbose)
    genome_file_paths = glob.glob(f"{get_tmp_dir(savedir)}/faa/*.faa")
    output_file = f"{get_tmp_dir(savedir)}/merged_genomes.faa"


    if os.path.exists(output_file):
        conditional_output(f"✅ [OK] Merged genomes file already exists at: {output_file}\n", "white", verbose)
        if logger:
            logger.info(f"Using existing merged genomes file: {output_file}")
        return True


    if not genome_file_paths:
        error_msg = f"No FAA files found in {get_tmp_dir(savedir)}/faa/"
        # Always show errors regardless of verbose setting
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        if logger:
            logger.error(error_msg)
        return False

    conditional_output("Merging individual genome files...", "yellow", verbose)
    try:
        with open(output_file, "w") as outfile:
            for each_file in genome_file_paths:
                with open(each_file) as infile:
                    for line in infile:
                        outfile.write(line)
        conditional_output(f"✅ [OK] Fasta files merged at: {output_file}\n", "white", verbose)
        if logger:
            logger.info(f"Successfully created merged genome file: {output_file}")
        return True
    except Exception as e:
        error_msg = f"Error merging genome files: {str(e)}"
        # Always show errors regardless of verbose setting
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        if logger:
            logger.error(error_msg)
        return False


def run_emapper(savedir: str, ncpus: int, resource_log_file: str, lowmem: bool = False, logger: Optional[logging.Logger] = None, verbose: bool = True) -> bool:
    """
    Run eggNOG-mapper on the merged genomes file.

    Parameters
    ----------
    savedir : str
        Directory where outputs will be saved
    ncpus : int
        Number of CPUs to use
    lowmem : bool, optional
        Whether to run with low memory settings, by default False
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    bool
        True if emapper ran successfully or outputs already exist, False otherwise
    """
    typer.secho("\n📊 Running emapper", fg="green")


    final_emapper_annotation_file = f"{savedir}/emapper_out.emapper.annotations"
    if os.path.exists(final_emapper_annotation_file):
        typer.secho(f"✅ [OK] Emapper annotations already exist at: {final_emapper_annotation_file}\n", fg="white")
        if logger:
            logger.info(f"Using existing emapper annotations: {final_emapper_annotation_file}")
        return True


    merged_genomes_file = f"{get_tmp_dir(savedir)}/merged_genomes.faa"
    output_folder_emapper = f"{get_tmp_dir(savedir)}/emapper_output"
    emapper_tmp_file = f"{output_folder_emapper}/emapper_out.emapper.annotations"


    if not os.path.exists(merged_genomes_file):
        error_msg = f"Merged genomes file not found at: {merged_genomes_file}"
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        if logger:
            logger.error(error_msg)
        return False


    if os.path.exists(emapper_tmp_file):
        typer.secho(f"✅ [OK] Emapper output already exists at: {emapper_tmp_file}\n", fg="white")
        if logger:
            logger.info(f"Using existing emapper output from temporary directory: {emapper_tmp_file}")

        try:
            shutil.copy(emapper_tmp_file, final_emapper_annotation_file)
            if logger:
                logger.info(f"Copied emapper annotations to: {final_emapper_annotation_file}")
        except Exception as e:
            if logger:
                logger.warning(f"Failed to copy emapper output to final location: {e}")
        return True


    if not os.path.exists(output_folder_emapper):
        os.makedirs(output_folder_emapper, exist_ok=True)


    try:
        cmd_emapper = [
            "stdbuf", "-o0", "-e0",  # Disable output buffering
            "emapper.py",
            "-i", merged_genomes_file,
            "--itype", "proteins",
            "--output_dir", output_folder_emapper,
            "--output", "emapper_out",
            "--cpu", str(ncpus)
        ]
        if not lowmem:
            cmd_emapper.append("--dbmem")

        returncode, stdout, stderr = run_subprocess_with_resource_monitoring(
            cmd_emapper,
            resource_log_file,
            logger,
            "Running eggNOG-mapper",
            verbose
        )

        if returncode != 0:
            error_msg = f"emapper failed with return code {returncode}"
            if stderr:
                error_msg += f": {stderr}"
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            if logger:
                logger.error(error_msg)
            return False


        if logger and stdout:
            logger.info(f"emapper stdout summary: {stdout[:500]}{'...' if len(stdout) > 500 else ''}")


        if not os.path.exists(emapper_tmp_file):
            error_msg = f"emapper did not generate expected output: {emapper_tmp_file}"
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            if logger:
                logger.error(error_msg)
            return False


        shutil.copy(emapper_tmp_file, final_emapper_annotation_file)

        typer.secho(f"✅ [OK] emapper output saved at: {output_folder_emapper}\n", fg="white")
        typer.secho(f"✅ [OK] emapper annotations copied to: {final_emapper_annotation_file}\n", fg="white")
        if logger:
            logger.info(f"Successfully ran emapper and saved annotations to: {final_emapper_annotation_file}")
        return True

    except Exception as e:
        error_msg = f"Error running emapper: {str(e)}"
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        if logger:
            logger.error(error_msg)
        return False


def remove_temp_files(savedir: str, logger: Optional[logging.Logger] = None) -> None:
    """
    Remove temporary files and directories.

    Parameters
    ----------
    savedir : str
        Directory where outputs are saved
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress
    """
    tmp_dir = get_tmp_dir(savedir)
    if os.path.exists(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
            typer.secho(f"✅ [OK] Temporary files removed from: {tmp_dir}", fg="white")
            if logger:
                logger.info(f"Removed temporary directory: {tmp_dir}")
        except Exception as e:
            typer.secho(f"⚠️ [WARNING] Failed to remove temporary files: {str(e)}", fg="yellow")
            if logger:
                logger.warning(f"Failed to remove temporary directory {tmp_dir}: {str(e)}")


def check_final_reports_exist(savedir: str, calculate_complementarity: int, logger: Optional[logging.Logger] = None) -> bool:
    """
    Check if all final output reports already exist.

    Parameters
    ----------
    savedir : str
        Directory where outputs are saved
    calculate_complementarity : int
        Maximum number of members for complementarity calculation
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    bool
        True if all required reports exist, False otherwise
    """

    ko_matrix_path = f"{savedir}/kos_matrix.csv"
    if not os.path.exists(ko_matrix_path):
        return False


    module_matrix_path = f"{savedir}/module_completeness.tsv"
    if not os.path.exists(module_matrix_path):
        return False


    if calculate_complementarity >= 2:
        for n_members in range(2, calculate_complementarity + 1):
            report_path = f"{savedir}/module_completeness_complementarity_{n_members}member.tsv"
            if not os.path.exists(report_path):
                return False


    if logger:
        logger.info(f"All required output files already exist in {savedir}")
    return True


def generate_complementarity_report(
    savedir: str,
    n_members: int = 2,
    logger: Optional[logging.Logger] = None
    ) -> None:
    """
    Generate a complementarity report showing modules completed only through genome partnerships.

    This function identifies metabolic modules that exhibit complementarity patterns:
    modules that are 100% complete when N genomes are combined but incomplete
    in all individual member genomes and smaller combinations. This reveals
    metabolic cooperation opportunities in microbial communities.

    The analysis works by:
    1. Reading module completeness data from the module matrix
    2. Parsing KPCT output to extract module metadata and KO mappings
    3. Optionally parsing emapper annotations for protein-level tracking
    4. Comparing completeness across individual genomes vs. combinations
    5. Identifying modules showing true complementarity patterns
    6. Generating a detailed report with module information and contributing proteins

    Parameters
    ----------
    savedir : str
        Directory containing module completeness matrix and where report will be saved
    n_members : int, optional
        Number of genomes in each combination to analyze, by default 2.
        Must be >= 2. Higher values find more complex complementarity patterns.
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress and statistics

    Output Files
    ------------
    module_completeness_complementarity_{n_members}member.tsv
        TSV file containing complementary modules with columns:
        - genome1, genome2, ..., genomeN: Individual genome IDs
        - genome1_completeness, ...: Completeness values for individuals
        - module_id: KEGG module identifier (e.g., M00001)
        - module_name: Human-readable module name
        - pathway_class: Metabolic pathway classification
        - matching_ko: KO terms involved in the module
        - genome1_proteins, ...: Contributing proteins from each genome (if available)

    Notes
    -----
    Requires pre-existing files:
    - module_completeness.tsv (from create_module_completeness_matrix)
    - KPCT output files for module metadata
    - Optionally: emapper annotations for protein tracking
    """


    output_file = f"{savedir}/module_completeness_complementarity_{n_members}member.tsv"


    if os.path.exists(output_file):
        if logger:
            logger.info(f"Complementarity report already exists at {output_file}")
        typer.secho(f"✅ [OK] Complementarity report already exists at: {output_file}", fg="white")
        return


    module_matrix_file = f"{savedir}/module_completeness.tsv"
    if not os.path.exists(module_matrix_file):
        error_msg = f"Module completeness matrix not found at: {module_matrix_file}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return


    emapper_file = None
    possible_emapper_files = [
        f"{savedir}/emapper_out.emapper.annotations",
        f"{get_tmp_dir(savedir)}/emapper_output/emapper_out.emapper.annotations",
        f"{savedir}/tmp/emapper_output/emapper_out.emapper.annotations",

    ]

    for possible_file in possible_emapper_files:
        if os.path.exists(possible_file):
            emapper_file = possible_file
            if logger:
                logger.info(f"Found emapper annotation file at: {emapper_file}")
            typer.secho(f"✅ [OK] Using emapper annotations from: {emapper_file}", fg="white")
            break

    if not emapper_file:
        if logger:
            logger.warning(f"Emapper annotation file not found in any of the expected locations. Will use placeholder protein IDs.")
        typer.secho(f"⚠️ [WARNING] Emapper annotation file not found. Will use placeholder protein IDs.", fg="yellow")


    kpct_output_file = None
    possible_kpct_files = [
        f"{savedir}/output_give_completeness_contigs.with_weights.tsv",
        f"{savedir}/output_give_completeness_pathways.with_weights.tsv",
        f"{savedir}/output_give_completeness_contigs.tsv",
        f"{savedir}/output_give_completeness_pathways.tsv"
    ]

    for file_path in possible_kpct_files:
        if os.path.exists(file_path):
            kpct_output_file = file_path
            break

    if not kpct_output_file:
        error_msg = "KPCT output file not found. Cannot extract module metadata."
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return

    try:

        if logger:
            logger.info(f"Reading module completeness matrix from {module_matrix_file}")

        module_df = pd.read_csv(module_matrix_file, sep='\t')


        if logger:
            logger.info(f"Reading KPCT output from {kpct_output_file}")

        kpct_df = pd.read_csv(kpct_output_file, sep='\t')


        genome_ko_proteins = {}
        if emapper_file:
            if logger:
                logger.info(f"Parsing protein information from {emapper_file}")
            genome_ko_proteins = parse_emapper_annotations(emapper_file, logger)


            if logger:
                logger.info(f"Found protein information for {len(genome_ko_proteins)} genomes")
                if genome_ko_proteins:
                    logger.debug(f"Genomes found in emapper annotations: {list(genome_ko_proteins.keys())}")


                taxon_oids = set(module_df[module_df['n_members'] == 1]['taxon_oid'])
                logger.debug(f"Found {len(taxon_oids)} unique genomes in module completeness matrix")


                if genome_ko_proteins:
                    emapper_genomes = set(genome_ko_proteins.keys())
                    missing_in_emapper = taxon_oids - emapper_genomes
                    if missing_in_emapper:
                        logger.warning(f"Found {len(missing_in_emapper)} genomes in module matrix but missing in emapper annotations")
                        if len(missing_in_emapper) < 10:
                            logger.warning(f"Missing genomes: {missing_in_emapper}")


        module_id_col = None
        module_name_col = None
        pathway_class_col = None
        contig_col = None
        matching_ko_col = None


        for col in ['module_accession', 'Module', 'module_id']:
            if col in kpct_df.columns:
                module_id_col = col
                break

        for col in ['pathway_name', 'Module_name', 'module_name']:
            if col in kpct_df.columns:
                module_name_col = col
                break

        for col in ['pathway_class', 'Module_class', 'class']:
            if col in kpct_df.columns:
                pathway_class_col = col
                break

        for col in ['contig', 'Contig', 'genome', 'Genome', 'taxon_oid']:
            if col in kpct_df.columns:
                contig_col = col
                break

        for col in ['matching_ko', 'module_ko', 'KOs']:
            if col in kpct_df.columns:
                matching_ko_col = col
                break

        if not module_id_col or not module_name_col:
            error_msg = f"Cannot identify required module columns in KPCT output: {kpct_output_file}"
            if logger:
                logger.error(error_msg)
                if module_id_col: logger.error(f"Found module_id_col: {module_id_col}")
                if module_name_col: logger.error(f"Found module_name_col: {module_name_col}")
                logger.error(f"Available columns: {kpct_df.columns.tolist()}")
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return

        if not contig_col or not matching_ko_col:
            if logger:
                logger.warning(f"Cannot identify contig or matching_ko columns in KPCT output.")
                logger.warning(f"Found contig_col: {contig_col}, matching_ko_col: {matching_ko_col}")
                logger.warning(f"Available columns: {kpct_df.columns.tolist()}")
            typer.secho(f"⚠️ [WARNING] Missing columns in KPCT output may affect mapping of KOs to combinations.", fg="yellow")


        module_metadata = {}
        module_to_kos = {}


        combo_module_kos = {}


        for _, row in kpct_df.iterrows():
            module_id = row[module_id_col]
            module_name = row[module_name_col]
            pathway_class = row[pathway_class_col] if pathway_class_col and pathway_class_col in row else ""


            if contig_col and matching_ko_col and not pd.isna(row[matching_ko_col]):
                taxon_oid = row[contig_col]
                ko_weights_str = row[matching_ko_col]
                if isinstance(ko_weights_str, str) and ko_weights_str:
                    combo_key = (taxon_oid, module_id)
                    combo_module_kos[combo_key] = ko_weights_str


            matching_kos = []
            for ko_col in ['matching_ko', 'module_ko', 'KOs']:
                if ko_col in kpct_df.columns and not pd.isna(row[ko_col]):
                    ko_str = row[ko_col]
                    if isinstance(ko_str, str) and ko_str:

                        for sep in [',', ' ', ';']:
                            if sep in ko_str:
                                matching_kos = [ko.strip() for ko in ko_str.split(sep) if ko.strip()]
                                break
                        if not matching_kos:
                            matching_kos = [ko_str.strip()]
                        break

            module_metadata[module_id] = {
                'module_name': module_name,
                'pathway_class': pathway_class,
            }

            if matching_kos:

                cleaned_kos = []
                for ko in matching_kos:

                    if '(' in ko:
                        ko = ko.split('(')[0].strip()
                    cleaned_kos.append(ko)

                module_to_kos[module_id] = cleaned_kos

        if not module_to_kos and logger:
            logger.warning("Could not extract KO lists for modules from KPCT output. Protein reporting will be limited.")
        else:
            if logger:
                logger.info(f"Extracted KO lists for {len(module_to_kos)} modules")

                sample_modules = list(module_to_kos.keys())[:3]
                for module in sample_modules:
                    logger.debug(f"Module {module} KOs: {module_to_kos[module]}")


                logger.info(f"Found KO weights for {len(combo_module_kos)} combination-module pairs")
                if combo_module_kos:

                    sample_keys = list(combo_module_kos.keys())[:3]
                    for key in sample_keys:
                        logger.debug(f"Combination-Module {key}: KOs with weights = {combo_module_kos[key]}")


        complementary_combinations = []


        single_genome_df = module_df[module_df['n_members'] == 1]
        n_member_df = module_df[module_df['n_members'] == n_members]

        if logger:
            logger.info(f"Found {len(single_genome_df)} single genomes and {len(n_member_df)} {n_members}-member combinations")


        module_cols = [col for col in module_df.columns if col.startswith('M') and col[1:].isdigit()]

        if logger:
            logger.info(f"Found {len(module_cols)} potential modules to evaluate")



        completeness_scale = 1.0
        for _, row in module_df.head().iterrows():
            for module_col in module_cols:
                if module_col in row and row[module_col] > 1.5:
                    completeness_scale = 100.0
                    break
            if completeness_scale == 100.0:
                break

        if logger:
            logger.info(f"Detected completeness scale: 0-{completeness_scale}")


        complete_in_combo_count = 0
        incomplete_in_individuals_count = 0
        incomplete_in_smaller_combos_count = 0


        for _, combo_row in n_member_df.iterrows():
            taxon_oids = combo_row['taxon_oid'].split('__')

            if len(taxon_oids) != n_members:
                continue


            individual_genomes = []
            for taxon_oid in taxon_oids:
                genome_row = single_genome_df[single_genome_df['taxon_oid'] == taxon_oid]
                if not genome_row.empty:
                    individual_genomes.append(genome_row.iloc[0])

            if len(individual_genomes) != n_members:
                continue

            for module_col in module_cols:

                if round(combo_row[module_col], 6) == completeness_scale:
                    complete_in_combo_count += 1


                    all_incomplete = all(round(genome[module_col], 6) < completeness_scale for genome in individual_genomes)

                    if all_incomplete:
                        incomplete_in_individuals_count += 1


                    smaller_combinations_incomplete = True
                    if n_members > 2:

                        for size in range(2, n_members):
                            for combo in itertools.combinations(taxon_oids, size):
                                combo_id = '__'.join(combo)
                                smaller_combo_row = module_df[module_df['taxon_oid'] == combo_id]
                                if not smaller_combo_row.empty and round(smaller_combo_row.iloc[0][module_col], 6) == completeness_scale:
                                    smaller_combinations_incomplete = False
                                    break
                            if not smaller_combinations_incomplete:
                                break

                    if all_incomplete and smaller_combinations_incomplete:
                        incomplete_in_smaller_combos_count += 1


                        module_kos = module_to_kos.get(module_col, [])


                        matching_ko_weights = ""
                        combo_id = combo_row['taxon_oid']
                        combo_module_key = (combo_id, module_col)
                        if combo_module_key in combo_module_kos:
                            matching_ko_weights = combo_module_kos[combo_module_key]


                        if not matching_ko_weights and contig_col and matching_ko_col:

                            module_rows = kpct_df[kpct_df[module_id_col] == module_col]

                            combo_row_kpct = module_rows[module_rows[contig_col] == combo_id]
                            if not combo_row_kpct.empty and matching_ko_col in combo_row_kpct.columns:
                                ko_value = combo_row_kpct.iloc[0][matching_ko_col]
                                if not pd.isna(ko_value):
                                    matching_ko_weights = ko_value


                        protein_contributions = {}
                        for taxon_oid in taxon_oids:
                            if genome_ko_proteins and module_kos:

                                if taxon_oid in genome_ko_proteins:

                                    genome_contributions = {}
                                    for ko in module_kos:
                                        if ko in genome_ko_proteins[taxon_oid]:
                                            genome_contributions[ko] = genome_ko_proteins[taxon_oid][ko]

                                    protein_contributions[taxon_oid] = genome_contributions


                                    if not genome_contributions and logger:
                                        logger.debug(f"Genome {taxon_oid} has protein data but no contributions to module {module_col}")
                                        logger.debug(f"Module KOs: {module_kos}")
                                        logger.debug(f"Genome KOs: {list(genome_ko_proteins[taxon_oid].keys())[:10]} (showing first 10)")
                                else:

                                    if logger:
                                        logger.debug(f"No protein data found for genome {taxon_oid}")
                                    protein_contributions[taxon_oid] = {}
                            else:

                                protein_contributions[taxon_oid] = {}


                        if logger:
                            logger.debug(f"Found complementary module: {module_col}")
                            logger.debug(f"  Combo completeness: {combo_row[module_col]}")
                            logger.debug(f"  Matching KOs with weights: {matching_ko_weights}")
                            for i, taxon_oid in enumerate(taxon_oids):
                                logger.debug(f"  {taxon_oid} completeness: {individual_genomes[i][module_col]}")
                                logger.debug(f"  {taxon_oid} contributing KOs: {list(protein_contributions.get(taxon_oid, {}).keys())}")


                        result = {
                            'module_id': module_col,
                            'module_name': module_metadata.get(module_col, {}).get('module_name', 'Unknown'),
                            'pathway_class': module_metadata.get(module_col, {}).get('pathway_class', ''),
                            'matching_ko': matching_ko_weights,
                        }


                        for i, taxon_oid in enumerate(taxon_oids):
                            result[f'taxon_oid_{i+1}'] = taxon_oid
                            genome_row = single_genome_df[single_genome_df['taxon_oid'] == taxon_oid]
                            result[f'completeness_taxon_oid_{i+1}'] = genome_row.iloc[0][module_col] if not genome_row.empty else 0


                            if taxon_oid in protein_contributions and protein_contributions[taxon_oid]:
                                result[f'proteins_taxon_oid_{i+1}'] = format_protein_contributions(taxon_oid, protein_contributions[taxon_oid])
                            else:

                                result[f'proteins_taxon_oid_{i+1}'] = f"No protein data available for {taxon_oid}"

                        complementary_combinations.append(result)


        if logger:
            logger.info(f"Modules exactly 100% complete in combinations: {complete_in_combo_count}")
            logger.info(f"Modules incomplete in all individuals: {incomplete_in_individuals_count}")
            if n_members > 2:
                logger.info(f"Modules incomplete in smaller combinations: {incomplete_in_smaller_combos_count}")
            logger.info(f"Total complementary modules found: {len(complementary_combinations)}")


        if complementary_combinations:
            report_df = pd.DataFrame(complementary_combinations)


            columns = []
            for i in range(1, n_members + 1):
                columns.append(f'taxon_oid_{i}')
            for i in range(1, n_members + 1):
                columns.append(f'completeness_taxon_oid_{i}')
            columns.extend(['module_id', 'module_name', 'pathway_class', 'matching_ko'])
            for i in range(1, n_members + 1):
                columns.append(f'proteins_taxon_oid_{i}')


            columns = [col for col in columns if col in report_df.columns]
            report_df = report_df[columns]


            report_df.to_csv(output_file, sep='\t', index=False)

            if logger:
                logger.info(f"Found {len(report_df)} complementary modules in {n_members}-member combinations")
                logger.info(f"Complementarity report saved to: {output_file}")

            typer.secho(f"✅ [OK] Found {len(report_df)} complementary modules in {n_members}-member combinations", fg="green")
            typer.secho(f"Complementarity report saved to: {output_file}", fg="white")
        else:
            if logger:
                logger.info(f"No modules found that are exactly 100% complete in {n_members}-member combinations but incomplete in all individuals")
                logger.info(f"Completeness scale detected: 0-{completeness_scale}")


            columns = []
            for i in range(1, n_members + 1):
                columns.append(f'taxon_oid_{i}')
            for i in range(1, n_members + 1):
                columns.append(f'completeness_taxon_oid_{i}')
            columns.extend(['module_id', 'module_name', 'pathway_class', 'matching_ko'])
            for i in range(1, n_members + 1):
                columns.append(f'proteins_taxon_oid_{i}')

            report_df = pd.DataFrame(columns=columns)
            report_df.to_csv(output_file, sep='\t', index=False)

            typer.secho(f"⚠️ [WARNING] No complementary modules found in {n_members}-member combinations", fg="yellow")
            typer.secho(f"Empty report saved to: {output_file}", fg="white")

    except Exception as e:
        error_msg = f"Error generating complementarity report: {str(e)}"
        if logger:
            logger.error(error_msg, exc_info=True)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")


def ko_matrix_to_kpct_format(kos_matrix: str, savedir: str, calculate_complementarity: int = 0, logger: Optional[logging.Logger] = None) -> str:
    """
    Convert KO matrix to KPCT format.

    Parameters
    ----------
    kos_matrix : str
        Path to the KO matrix CSV/TSV file
    savedir : str
        Directory to save the KPCT format file
    calculate_complementarity : int, optional
        Number of members for complementarity calculation, by default 0
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    str
        Path to the generated KPCT format file
    """
    if logger:
        logger.info(f"Converting KO matrix to KPCT format: {kos_matrix}")

    output_path = os.path.join(savedir, "ko_file_for_kpct.txt")

    try:

        initial_delimiter = ',' if kos_matrix.lower().endswith('.csv') else '\t'

        typer.secho(f"Reading KO matrix file: {kos_matrix}", fg="yellow")
        if logger:
            logger.info(f"Reading KO matrix file with delimiter '{initial_delimiter}': {kos_matrix}")


        try:
            ko_df = pd.read_csv(kos_matrix, sep=initial_delimiter)
        except pd.errors.EmptyDataError as e_empty:
            if kos_matrix.lower().endswith('.tsv') and initial_delimiter == '\t':
                if logger:
                    logger.info(f"File is .tsv and tab-delimited read failed. Attempting fallback to comma delimiter for {kos_matrix}.")
                try:
                    ko_df = pd.read_csv(kos_matrix, sep=',')
                    if logger:
                        logger.info(f"Successfully read {kos_matrix} with fallback comma delimiter.")
                except Exception as e_fallback:
                    if logger:
                        logger.error(f"Fallback to comma delimiter for {kos_matrix} also failed: {e_fallback}")

                    try:
                        file_size = os.path.getsize(kos_matrix)
                        if file_size == 0:
                            if logger:
                                logger.error(f"Confirmed: File {kos_matrix} is 0 bytes.")
                        else:
                            if logger:
                                logger.error(f"File {kos_matrix} is {file_size} bytes but still failed to parse with fallback.")
                    except OSError as e_os:
                        if logger:
                            logger.error(f"Could not get size of file {kos_matrix}: {e_os}")
                    raise e_empty
            else:
                try:
                    if os.path.getsize(kos_matrix) == 0 and logger:
                        logger.error(f"File {kos_matrix} is 0 bytes.")
                except OSError: pass
                raise e_empty
        except Exception as e_general:
            if logger:
                logger.error(f"General error reading {kos_matrix} with delimiter '{initial_delimiter}': {e_general}")
            raise e_general

        if ko_df is None:
             final_error_msg = f"Failed to load DataFrame from {kos_matrix} after all attempts; ko_df is None."
             if logger:
                 logger.error(final_error_msg)

             raise pd.errors.EmptyDataError(final_error_msg)



        if 'taxon_oid' not in ko_df.columns:
            msg = "Invalid KO matrix format: missing 'taxon_oid' column"
            if logger:
                logger.error(msg)
            typer.secho(f"❌ [ERROR] {msg}", fg="red")
            exit(1)


        taxon_oids = ko_df['taxon_oid'].tolist()
        ko_columns = [col for col in ko_df.columns if col != 'taxon_oid']


        with open(output_path, 'w') as f:

            for idx, row in ko_df.iterrows():
                taxon_id = row['taxon_oid']
                present_kos = []

                for ko_id in ko_columns:
                    if row[ko_id] > 0:
                        present_kos.append(ko_id)

                if present_kos:
                    f.write(f"{taxon_id}\t" + "\t".join(present_kos) + "\n")


            if calculate_complementarity >= 2:

                combinations = []
                for n in range(2, calculate_complementarity + 1):
                    combinations.extend(list(itertools.combinations(taxon_oids, n)))


                for combo in combinations:
                    combo_taxon_id = "__".join(combo)
                    combo_kos = set()


                    for taxon_id in combo:
                        taxon_row = ko_df[ko_df['taxon_oid'] == taxon_id].iloc[0]
                        for ko_id in ko_columns:
                            if taxon_row[ko_id] > 0:
                                combo_kos.add(ko_id)

                    if combo_kos:
                        f.write(f"{combo_taxon_id}\t" + "\t".join(sorted(combo_kos)) + "\n")

        if logger:
            logger.info(f"KO matrix converted to KPCT format: {output_path}")
        typer.secho(f"✅ [OK] KO matrix converted to KPCT format: {output_path}", fg="white")
        return output_path

    except Exception as e:
        error_msg = f"Error converting KO matrix to KPCT format: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        raise


def get_ko_protein_mappings_from_kpct_input(kpct_input_file: str, logger: Optional[logging.Logger] = None) -> Dict[str, Dict[str, str]]:
    """
    Parse the KPCT input file to extract KO to protein mappings for each genome.

    Parameters
    ----------
    kpct_input_file : str
        Path to the KPCT input file which contains genome/KO mappings
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    Dict[str, Dict[str, str]]
        A dictionary mapping genome IDs to dictionaries of KO IDs to protein IDs
        {genome_id: {ko_id: protein_id, ...}, ...}
    """
    ko_protein_mappings = {}

    try:
        with open(kpct_input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    continue

                genome_id = parts[0]
                ko_ids = parts[1:]



                ko_protein_map = {}
                for ko_id in ko_ids:
                    protein_id = f"{genome_id}|{ko_id}"
                    ko_protein_map[ko_id] = protein_id

                ko_protein_mappings[genome_id] = ko_protein_map

        if logger:
            logger.info(f"Extracted KO-protein mappings for {len(ko_protein_mappings)} genomes from KPCT input file")

        return ko_protein_mappings
    except Exception as e:
        if logger:
            logger.error(f"Error extracting KO-protein mappings from KPCT input file: {str(e)}")
        return {}


def create_module_completeness_matrix(savedir: str, kpct_outprefix: str, logger: Optional[logging.Logger] = None) -> None:
    """
    Create a module completeness matrix from the KPCT output.

    Parameters
    ----------
    savedir : str
        Directory where KPCT output is saved and where to write the matrix
    kpct_outprefix : str
        Prefix for KPCT output files
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress
    """
    output_file = os.path.join(savedir, "module_completeness.tsv")
    kpct_output_file = os.path.join(savedir, f"{kpct_outprefix}_contigs.with_weights.tsv")

    # Try to find the KPCT input file to get the complete list of genomes
    kpct_input_file = os.path.join(savedir, "ko_file_for_kpct.txt")

    if not os.path.exists(kpct_output_file):
        alternatives = [
            os.path.join(savedir, f"{kpct_outprefix}_pathways.with_weights.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_contigs.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_pathways.tsv")
        ]
        for alt_file in alternatives:
            if os.path.exists(alt_file):
                kpct_output_file = alt_file
                break

    if not os.path.exists(kpct_output_file):
        error_msg = f"KPCT output file not found: tried {kpct_outprefix}_contigs.with_weights.tsv and alternatives"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return

    try:
        # Get the complete list of genomes from the original KPCT input
        all_genomes = set()
        if os.path.exists(kpct_input_file):
            with open(kpct_input_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        genome_id = line.split('\t')[0]
                        all_genomes.add(genome_id)
            if logger:
                logger.info(f"Found {len(all_genomes)} total genomes in KPCT input file")

        if logger:
            logger.info(f"Reading KPCT output from: {kpct_output_file}")

        kpct_df = pd.read_csv(kpct_output_file, sep='\t')

        # Identify contig/genome column
        contig_column = None
        possible_contig_columns = ['contig', 'Contig', 'genome', 'Genome', 'taxon_oid']

        for col in possible_contig_columns:
            if col in kpct_df.columns:
                contig_column = col
                break

        if not contig_column:
            # Fall back to first column
            contig_column = kpct_df.columns[0]
            if logger:
                logger.warning(f"No standard contig/genome column found. Using first column: '{contig_column}'")

        # Identify module columns and completeness column
        module_cols = []
        completeness_col = None

        for col in kpct_df.columns:
            if col.startswith('M') and col[1:].isdigit():
                module_cols.append(col)
            elif col.lower() == 'completeness':
                completeness_col = col

        # Handle row-based format (pivot if needed)
        if not module_cols and 'module_accession' in kpct_df.columns and completeness_col:
            if logger:
                logger.info("Detected row-based format. Pivoting data to create module matrix.")

            # Get all unique modules before pivoting
            all_modules = sorted(kpct_df['module_accession'].unique())

            # Pivot the data
            pivot_df = kpct_df.pivot_table(
                index=contig_column,
                columns='module_accession',
                values=completeness_col,
                fill_value=0
            )

            # Reset index to make contig_column a regular column
            pivot_df = pivot_df.reset_index()

            # If we have a list of all genomes, ensure they're all included
            if all_genomes:
                # Get genomes that appear in the KPCT output
                genomes_in_output = set(pivot_df[contig_column].astype(str))

                # Find missing single genomes (not combinations)
                missing_single_genomes = set()
                for genome in all_genomes:
                    if genome not in genomes_in_output and '__' not in genome:
                        missing_single_genomes.add(genome)

                if missing_single_genomes:
                    if logger:
                        logger.info(f"Adding {len(missing_single_genomes)} missing genomes with zero completeness")

                    # Create rows for missing genomes with zero completeness
                    missing_rows = []
                    for genome in missing_single_genomes:
                        row = {contig_column: genome}
                        for module in all_modules:
                            row[module] = 0
                        missing_rows.append(row)

                    # Add missing rows to the pivot dataframe
                    if missing_rows:
                        missing_df = pd.DataFrame(missing_rows)
                        pivot_df = pd.concat([pivot_df, missing_df], ignore_index=True)

            kpct_df = pivot_df
            module_cols = [col for col in kpct_df.columns if col != contig_column]

        if not module_cols:
            error_msg = "Could not identify module columns in the KPCT output"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return

        # Build the result data
        result_data = []

        for _, row in kpct_df.iterrows():
            taxon_oid = row[contig_column]

            # Determine number of members
            n_members = 1
            if "__" in str(taxon_oid):
                n_members = str(taxon_oid).count("__") + 1

            # Create result row
            result_row = {'n_members': n_members, 'taxon_oid': taxon_oid}

            # Add module completeness values
            for module in module_cols:
                result_row[module] = row[module]

            result_data.append(result_row)

        # Create final dataframe
        result_df = pd.DataFrame(result_data)

        # Sort columns: n_members, taxon_oid, then modules alphabetically
        module_cols = sorted([col for col in result_df.columns if col not in ['n_members', 'taxon_oid']])
        result_df = result_df[['n_members', 'taxon_oid'] + module_cols]

        # Sort by n_members and taxon_oid
        result_df = result_df.sort_values(by=['n_members', 'taxon_oid'])

        # Save the matrix
        result_df.to_csv(output_file, sep='\t', index=False)

        if logger:
            single_genomes = len(result_df[result_df['n_members'] == 1])
            total_genomes = len(result_df)
            logger.info(f"Module completeness matrix saved to: {output_file}")
            logger.info(f"Matrix contains {single_genomes} single genomes out of {total_genomes} total entries")
            if all_genomes:
                logger.info(f"Expected {len(all_genomes)} single genomes from KPCT input")
        typer.secho(f"✅ [OK] Module completeness matrix saved to: {output_file}", fg="white")

    except Exception as e:
        error_msg = f"Error creating module completeness matrix: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")

        if logger:
            logger.error(f"Error details: {e}", exc_info=True)


def create_ko_matrix_from_emapper_annotation(emapper_file_path: str, output_file_path: str, logger: Optional[logging.Logger] = None) -> None:
    """
    Create a KO matrix from an eggNOG-mapper annotation file.

    This function processes eggNOG-mapper output to create a genome × KO count matrix
    suitable for downstream module completeness analysis. The process includes:
    1. Reading eggNOG-mapper annotations (skipping header/footer lines)
    2. Extracting genome IDs from protein headers (format: genome_id|protein_id)
    3. Parsing KO annotations and cleaning KO identifiers
    4. Counting KO occurrences per genome
    5. Creating a pivot table with genomes as rows and KOs as columns
    6. Saving as CSV file for downstream processing

    Parameters
    ----------
    emapper_file_path : str
        Path to the emapper.annotations file produced by eggNOG-mapper.
        Expected to have tab-separated columns including "#query" and "KEGG_ko"
    output_file_path : str
        Path where the KO matrix CSV will be saved.
        Output format: rows=genomes, columns=KO terms, values=protein counts
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress and statistics

    Output Format
    -------------
    CSV file with:
    - Column 1: 'taxon_oid' containing genome identifiers
    - Remaining columns: KO identifiers (K00001, K00002, etc.)
    - Values: Integer counts of proteins for each KO in each genome
    - Missing KOs have value 0

    Notes
    -----
    - Assumes protein headers follow format: genome_id|protein_id
    - Handles comma-separated KO lists in eggNOG-mapper output
    - Removes 'ko:' prefixes and weight annotations like '(0.5)'
    - Skips rows with missing or '-' KO annotations
    """
    typer.secho("\n📊 Creating KO matrix from eggNOG-mapper annotations", fg="green")


    if os.path.exists(output_file_path):
        typer.secho(f"✅ [OK] KO matrix already exists at: {output_file_path}", fg="white")
        if logger:
            logger.info(f"KO matrix already exists at: {output_file_path}")
        return

    if not os.path.exists(emapper_file_path):
        error_msg = f"eMapper annotation file not found at {emapper_file_path}. Cannot proceed."
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        exit(1)

    try:
        typer.secho("Processing eggNOG-mapper annotations and extracting KO terms...", fg="yellow")

        if logger:
            logger.info(f"Reading eggNOG-mapper annotations from: {emapper_file_path}")



        emapper_df = pd.read_csv(
            emapper_file_path,
            sep="\t",
            skiprows=4,
            skipfooter=3,
            engine='python',
            usecols=["#query", "KEGG_ko"]
        )


        emapper_df = emapper_df[emapper_df["KEGG_ko"] != "-"]
        emapper_df = emapper_df.dropna(subset=["KEGG_ko", "#query"])
        emapper_df = emapper_df.reset_index(drop=True)

        if logger:
            logger.info(f"Loaded {len(emapper_df)} rows with KO annotations")


        emapper_df["genome_id"] = emapper_df["#query"].str.split('|').str[0]
        taxa_list = list(emapper_df["genome_id"].unique())

        if logger:
            logger.info(f"Found {len(taxa_list)} unique genomes in annotations")


        kos_data_for_matrix_df = []


        for each_taxon in taxa_list:
            genome_specific_df = emapper_df[emapper_df["genome_id"] == each_taxon]


            kos_in_genome = []

            for _, row in genome_specific_df.iterrows():
                kegg_ko_str = row["KEGG_ko"]
                if pd.isna(kegg_ko_str) or not kegg_ko_str:
                    continue


                row_list = kegg_ko_str.split(",")
                for each_nested_ko_raw in row_list:
                    ko_id = each_nested_ko_raw.strip().removeprefix("ko:").strip()
                    if ko_id:
                        kos_in_genome.append(ko_id)


            tmp_df_matrix = pd.DataFrame(kos_in_genome, columns=["ko"])
            tmp_df_matrix["taxon_oid"] = each_taxon
            kos_data_for_matrix_df.append(tmp_df_matrix)

        if not kos_data_for_matrix_df:
            error_msg = "No KO data found in the eMapper annotations file"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return


        typer.secho("Creating KO count matrix (kos_matrix.csv)...", fg="yellow")
        kos_count_df = pd.concat(kos_data_for_matrix_df)


        kos_count_df = kos_count_df.value_counts().rename_axis(["ko", "taxon_oid"]).reset_index(name="counts")
        kos_count_df['counts'] = kos_count_df['counts'].astype(int)


        kos_count_df = kos_count_df.pivot_table(
            index=["taxon_oid"],
            columns="ko",
            values="counts",
            fill_value=0
        ).reset_index(drop=False)


        for col in kos_count_df.columns:
            if col != 'taxon_oid':
                kos_count_df[col] = kos_count_df[col].astype(int)


        kos_count_df.to_csv(output_file_path, index=False)

        if logger:
            logger.info(f"Created KO matrix with {len(kos_count_df)} genomes and {len(kos_count_df.columns)-1} KOs")
            logger.info(f"KO matrix saved to: {output_file_path}")

        typer.secho(f"✅ [OK] KO matrix created and saved to: {output_file_path}", fg="white")

    except Exception as e:
        error_msg = f"Error creating KO matrix: {str(e)}"
        if logger:
            logger.error(error_msg, exc_info=True)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        exit(1)


def check_kpct_installed(logger: Optional[logging.Logger] = None) -> bool:
    """
    Check if the KPCT give_completeness tool is installed.

    Parameters
    ----------
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    bool
        True if KPCT is installed, False otherwise
    """
    try:

        give_completeness_check = subprocess.run(
            ["which", "give_completeness"],
            capture_output=True,
            text=True
        )

        if give_completeness_check.returncode != 0:
            error_msg = "KPCT 'give_completeness' tool not found in PATH. Please install it via pip: pip install kegg-pathways-completeness"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return False
        return True
    except Exception as e:
        error_msg = f"Error checking for KPCT installation: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return False


def check_chunk_outputs_exist(chunk_dirs: List[str], kpct_outprefix: str, logger: Optional[logging.Logger] = None) -> Tuple[bool, List[str]]:
    """
    Check if chunk outputs already exist and which ones are complete.

    Parameters
    ----------
    chunk_dirs : List[str]
        List of chunk directories to check
    kpct_outprefix : str
        Prefix for KPCT output files
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    Tuple[bool, List[str]]
        (All chunks exist, List of existing chunk directories)
    """
    existing_chunk_dirs = []

    for chunk_dir in chunk_dirs:
        if os.path.exists(chunk_dir):

            expected_files = [
                f"*_contigs.with_weights.tsv",
                f"*_pathways.with_weights.tsv"
            ]

            files_exist = True
            for pattern in expected_files:
                matching_files = glob.glob(os.path.join(chunk_dir, pattern))
                if not matching_files:
                    files_exist = False
                    break

            if files_exist:
                existing_chunk_dirs.append(chunk_dir)
                if logger:
                    logger.info(f"Found existing chunk output: {chunk_dir}")

    all_exist = len(existing_chunk_dirs) == len(chunk_dirs)

    if logger:
        logger.info(f"Chunk status: {len(existing_chunk_dirs)}/{len(chunk_dirs)} chunks already processed")

    return all_exist, existing_chunk_dirs


def chunk_kpct_input_file(kpct_input_file: str, savedir: str, n_chunks: int, logger: Optional[logging.Logger] = None) -> List[str]:
    """
    Split the KPCT input file into chunks for parallel processing.

    Parameters
    ----------
    kpct_input_file : str
        Path to the KPCT input file
    savedir : str
        Directory to save chunk files
    n_chunks : int
        Number of chunks to create
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    List[str]
        List of paths to chunk files
    """
    chunks_dir = os.path.join(get_tmp_dir(savedir), "kpct_chunks")
    os.makedirs(chunks_dir, exist_ok=True)

    if logger:
        logger.info(f"Chunking KPCT input file into {n_chunks} chunks")

    try:

        with open(kpct_input_file, 'r') as f:
            lines = f.readlines()

        if not lines:
            error_msg = f"KPCT input file is empty: {kpct_input_file}"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return []

        # Calculate lines per chunk using ceiling division to ensure we create exactly n_chunks
        # when there are enough lines, or fewer chunks only if we have fewer lines than chunks
        if len(lines) < n_chunks:
            # If we have fewer lines than requested chunks, create one chunk per line
            n_chunks = len(lines)
            if logger:
                logger.warning(f"Only {len(lines)} lines available, reducing chunks to {n_chunks}")

        lines_per_chunk = len(lines) // n_chunks
        remainder = len(lines) % n_chunks

        chunk_files = []

        # Create chunks with more even distribution
        current_idx = 0
        for i in range(n_chunks):
            # Some chunks get an extra line if there's a remainder
            chunk_size = lines_per_chunk + (1 if i < remainder else 0)

            # Skip if we've already processed all lines
            if current_idx >= len(lines):
                break

            chunk_lines = lines[current_idx:current_idx + chunk_size]

            if chunk_lines:
                chunk_file = os.path.join(chunks_dir, f"chunk_{i:03d}.txt")
                with open(chunk_file, 'w') as f:
                    f.writelines(chunk_lines)
                chunk_files.append(chunk_file)

            current_idx += chunk_size

        if logger:
            # Calculate actual distribution for better logging
            total_lines = len(lines)
            chunks_created = len(chunk_files)
            base_lines = lines_per_chunk
            extra_lines_chunks = remainder

            if extra_lines_chunks > 0:
                logger.info(f"Created {chunks_created} chunks: {extra_lines_chunks} chunks with {base_lines + 1} lines, {chunks_created - extra_lines_chunks} chunks with {base_lines} lines (total: {total_lines} KPCT input lines)")
            else:
                logger.info(f"Created {chunks_created} chunks with {base_lines} lines each (total: {total_lines} KPCT input lines)")

        return chunk_files

    except Exception as e:
        error_msg = f"Error chunking KPCT input file: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return []


def run_kpct_on_chunk(chunk_file: str, chunk_savedir: str, kpct_outprefix: str, chunk_id: str) -> Tuple[bool, str, str]:
    """
    Run KPCT on a single chunk file.

    Parameters
    ----------
    chunk_file : str
        Path to the chunk file
    chunk_savedir : str
        Directory to save chunk outputs
    kpct_outprefix : str
        Prefix for KPCT output files
    chunk_id : str
        Identifier for this chunk

    Returns
    -------
    Tuple[bool, str, str]
        Success status, chunk output directory, and error message if any
    """
    try:

        os.makedirs(chunk_savedir, exist_ok=True)


        chunk_prefix = f"{kpct_outprefix}_chunk_{chunk_id}"
        kpct_cmd = [
            "give_completeness",
            "--input", chunk_file,
            "--outprefix", chunk_prefix,
            "--include-weights",
            "--add-per-contig",
            "--outdir", chunk_savedir
        ]


        process = subprocess.run(
            kpct_cmd,
            capture_output=True,
            text=True

        )

        if process.returncode != 0:
            error_msg = f"KPCT failed for chunk {chunk_id} with return code {process.returncode}"
            if process.stderr:
                error_msg += f": {process.stderr}"
            return False, chunk_savedir, error_msg


        expected_files = [
            os.path.join(chunk_savedir, f"{chunk_prefix}_contigs.with_weights.tsv"),
            os.path.join(chunk_savedir, f"{chunk_prefix}_pathways.with_weights.tsv")
        ]

        files_created = [f for f in expected_files if os.path.exists(f)]
        if not files_created:
            return False, chunk_savedir, f"No output files created for chunk {chunk_id}"

        return True, chunk_savedir, ""

    except subprocess.TimeoutExpired:
        return False, chunk_savedir, f"KPCT timed out for chunk {chunk_id}"
    except Exception as e:
        return False, chunk_savedir, f"Error processing chunk {chunk_id}: {str(e)}"


def concatenate_kpct_outputs(chunk_dirs: List[str], savedir: str, kpct_outprefix: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    Concatenate KPCT outputs from all chunks into final output files.

    Parameters
    ----------
    chunk_dirs : List[str]
        List of directories containing chunk outputs
    savedir : str
        Directory to save final concatenated outputs
    kpct_outprefix : str
        Prefix for final KPCT output files
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    bool
        True if concatenation was successful, False otherwise
    """
    try:
        output_types = ["contigs.with_weights", "pathways.with_weights"]

        for output_type in output_types:
            final_output_file = os.path.join(savedir, f"{kpct_outprefix}_{output_type}.tsv")

            if logger:
                logger.info(f"Concatenating {output_type} outputs to {final_output_file}")

            chunk_files = []
            header_written = False


            for chunk_dir in chunk_dirs:
                pattern = os.path.join(chunk_dir, f"*_chunk_*_{output_type}.tsv")
                chunk_files.extend(glob.glob(pattern))

            if not chunk_files:
                if logger:
                    logger.warning(f"No chunk files found for {output_type}")
                continue


            chunk_files.sort()


            with open(final_output_file, 'w') as outfile:
                for chunk_file in chunk_files:
                    if os.path.exists(chunk_file):
                        with open(chunk_file, 'r') as infile:
                            lines = infile.readlines()

                            if lines:

                                if not header_written:
                                    outfile.writelines(lines)
                                    header_written = True
                                else:

                                    outfile.writelines(lines[1:])

            if logger:
                logger.info(f"Successfully created {final_output_file}")

        return True

    except Exception as e:
        error_msg = f"Error concatenating KPCT outputs: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return False


def run_kpct_parallel(kpct_input_file: str, savedir: str, kpct_outprefix: str, ncpus: int, logger: Optional[logging.Logger] = None) -> bool:
    """
    Run KPCT in parallel by chunking the input file and processing chunks concurrently.

    This function implements parallel processing for KPCT to improve performance on large datasets.
    The workflow includes:
    1. Check if final outputs already exist (skip if so)
    2. Split input file into chunks equal to number of CPU cores
    3. Check which chunks have already been processed (checkpointing)
    4. Process remaining chunks in parallel using ProcessPoolExecutor
    5. Concatenate all chunk outputs into final result files
    6. Validate that all expected outputs were created

    The parallel approach can significantly reduce processing time for large datasets
    while providing fault tolerance through checkpointing of individual chunks.

    Parameters
    ----------
    kpct_input_file : str
        Path to the KPCT input file containing genome/KO mappings
    savedir : str
        Directory to save outputs and temporary chunk files
    kpct_outprefix : str
        Prefix for KPCT output files (e.g., "output_give_completeness")
    ncpus : int
        Number of CPU cores to use (determines number of chunks created)
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress and debugging information

    Returns
    -------
    bool
        True if KPCT parallel processing completed successfully, False otherwise

    Notes
    -----
    - Creates temporary directories: tmp/kpct_chunks/ and tmp/kpct_chunk_outputs/
    - Supports checkpointing: existing chunk outputs are reused
    - Falls back gracefully if chunks cannot be processed
    - Concatenates results while preserving headers correctly
    """
    try:

        final_outputs = [
            os.path.join(savedir, f"{kpct_outprefix}_contigs.with_weights.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_pathways.with_weights.tsv")
        ]

        if all(os.path.exists(f) for f in final_outputs):
            if logger:
                logger.info("KPCT output files already exist, skipping parallel processing")
            typer.secho("✅ [OK] KPCT output files already exist", fg="white")
            return True


        n_chunks = ncpus

        if logger:
            logger.info(f"Running KPCT in parallel with up to {n_chunks} chunks using {ncpus} CPU cores")

        typer.secho(f"🔧 Running KPCT in parallel with up to {n_chunks} chunks", fg="yellow")


        chunks_base_dir = os.path.join(get_tmp_dir(savedir), "kpct_chunk_outputs")
        os.makedirs(chunks_base_dir, exist_ok=True)


        expected_chunk_dirs = []
        for i in range(n_chunks):
            chunk_id = f"{i:03d}"
            chunk_savedir = os.path.join(chunks_base_dir, f"chunk_{chunk_id}")
            expected_chunk_dirs.append(chunk_savedir)


        all_chunks_exist, existing_chunk_dirs = check_chunk_outputs_exist(expected_chunk_dirs, kpct_outprefix, logger)

        if all_chunks_exist:
            if logger:
                logger.info("All chunk outputs already exist, proceeding to concatenation")
            typer.secho("✅ [OK] All chunks already processed, concatenating results", fg="white")


            concatenation_success = concatenate_kpct_outputs(existing_chunk_dirs, savedir, kpct_outprefix, logger)
            return concatenation_success


        chunks_dir = os.path.join(get_tmp_dir(savedir), "kpct_chunks")
        chunk_files = []


        for i in range(n_chunks):
            chunk_file = os.path.join(chunks_dir, f"chunk_{i:03d}.txt")
            if os.path.exists(chunk_file):
                chunk_files.append(chunk_file)


        if len(chunk_files) != n_chunks:
            if logger:
                logger.info(f"Creating chunk files (found {len(chunk_files)}/{n_chunks} existing)")
            chunk_files = chunk_kpct_input_file(kpct_input_file, savedir, n_chunks, logger)
        else:
            if logger:
                logger.info(f"Using existing {len(chunk_files)} chunk files")

        if not chunk_files:
            return False

        # Update n_chunks and expected_chunk_dirs to match actual chunks created
        actual_n_chunks = len(chunk_files)
        if actual_n_chunks != n_chunks:
            if logger:
                logger.info(f"Adjusting expected chunks from {n_chunks} to {actual_n_chunks} based on actual data")
            n_chunks = actual_n_chunks
            # Recreate expected_chunk_dirs with the correct number
            expected_chunk_dirs = []
            for i in range(n_chunks):
                chunk_id = f"{i:03d}"
                chunk_savedir = os.path.join(chunks_base_dir, f"chunk_{chunk_id}")
                expected_chunk_dirs.append(chunk_savedir)


        chunks_to_process = []
        chunk_dirs_map = {}

        for i, chunk_file in enumerate(chunk_files):
            chunk_id = f"{i:03d}"
            chunk_savedir = os.path.join(chunks_base_dir, f"chunk_{chunk_id}")
            chunk_dirs_map[chunk_id] = chunk_savedir


            if chunk_savedir not in existing_chunk_dirs:
                chunks_to_process.append((chunk_file, chunk_savedir, chunk_id))

        if not chunks_to_process:
            if logger:
                logger.info("All chunks already processed, proceeding to concatenation")
            typer.secho("✅ [OK] All chunks already processed, concatenating results", fg="white")
        else:
            if logger:
                logger.info(f"Processing {len(chunks_to_process)} remaining chunks")
            typer.secho(f"🔧 Processing {len(chunks_to_process)} remaining chunks", fg="yellow")


            failed_chunks = []

            with ProcessPoolExecutor(max_workers=ncpus) as executor:

                future_to_chunk = {}

                for chunk_file, chunk_savedir, chunk_id in chunks_to_process:
                    future = executor.submit(
                        run_kpct_on_chunk,
                        chunk_file,
                        chunk_savedir,
                        kpct_outprefix,
                        chunk_id
                    )
                    future_to_chunk[future] = (chunk_id, chunk_savedir)


                for future in as_completed(future_to_chunk):
                    chunk_id, chunk_savedir = future_to_chunk[future]

                    try:
                        success, output_dir, error_msg = future.result()

                        if success:
                            if logger:
                                logger.info(f"Successfully processed chunk {chunk_id}")
                        else:
                            failed_chunks.append((chunk_id, error_msg))
                            if logger:
                                logger.error(f"Failed to process chunk {chunk_id}: {error_msg}")

                    except Exception as e:
                        failed_chunks.append((chunk_id, str(e)))
                        if logger:
                            logger.error(f"Exception processing chunk {chunk_id}: {str(e)}")


            if failed_chunks:
                error_msg = f"Failed to process {len(failed_chunks)} chunks: {failed_chunks}"
                if logger:
                    logger.error(error_msg)
                typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
                return False


        all_chunk_dirs = list(expected_chunk_dirs)


        all_chunks_exist, final_chunk_dirs = check_chunk_outputs_exist(all_chunk_dirs, kpct_outprefix, logger)

        if not all_chunks_exist:
            error_msg = f"Not all chunks were processed successfully. Expected {len(all_chunk_dirs)}, got {len(final_chunk_dirs)}"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return False


        if logger:
            logger.info(f"Concatenating outputs from {len(final_chunk_dirs)} chunks")

        concatenation_success = concatenate_kpct_outputs(final_chunk_dirs, savedir, kpct_outprefix, logger)

        if not concatenation_success:
            return False


        missing_outputs = [f for f in final_outputs if not os.path.exists(f)]
        if missing_outputs:
            error_msg = f"Failed to create final output files: {missing_outputs}"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return False

        if logger:
            logger.info("Successfully completed parallel KPCT processing")
        typer.secho("✅ [OK] KPCT parallel processing completed successfully", fg="green")

        return True

    except Exception as e:
        error_msg = f"Error in parallel KPCT processing: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return False


def run_kpct(kpct_input_file: str, savedir: str, kpct_outprefix: str, resource_log_file: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    Run the KPCT give_completeness tool (sequential version).
    This function is kept as a fallback in case parallel processing fails.

    Parameters
    ----------
    kpct_input_file : str
        Path to the KPCT input file
    savedir : str
        Directory to save outputs
    kpct_outprefix : str
        Prefix for KPCT output files
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    bool
        True if KPCT ran successfully, False otherwise
    """
    try:

        kpct_cmd = [
            "give_completeness",
            "--input", kpct_input_file,
            "--outprefix", kpct_outprefix,
            "--include-weights",
            "--add-per-contig",
            "--outdir", savedir
        ]


        returncode, stdout, stderr = run_subprocess_with_resource_monitoring(
            kpct_cmd,
            resource_log_file,
            logger,
            "Running KPCT give_completeness tool (sequential)",
            True
        )

        if returncode != 0:
            error_msg = f"KPCT tool failed with return code {returncode}"
            if stderr:
                error_msg += f": {stderr}"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return False


        if logger and stdout:
            logger.info(f"KPCT stdout: {stdout}")


        possible_kpct_files = [
            os.path.join(savedir, f"{kpct_outprefix}_contigs.with_weights.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_pathways.with_weights.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_contigs.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_pathways.tsv")
        ]
        kpct_file_exists = any(os.path.exists(f) for f in possible_kpct_files)

        if not kpct_file_exists:
            error_msg = f"KPCT did not generate any output files with prefix '{kpct_outprefix}'"
            if logger:
                logger.error(error_msg)
            typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
            return False


        created_files = [f for f in possible_kpct_files if os.path.exists(f)]
        if logger:
            logger.info(f"KPCT successfully created output files: {created_files}")
        typer.secho(f"✅ [OK] KPCT completed successfully. Created files: {[os.path.basename(f) for f in created_files]}", fg="green")

        return True

    except Exception as e:
        error_msg = f"Error running KPCT: {str(e)}"
        if logger:
            logger.error(error_msg)
        typer.secho(f"❌ [ERROR] {error_msg}", fg="red")
        return False


def run_kpct_with_fallback(kpct_input_file: str, savedir: str, kpct_outprefix: str, ncpus: int, resource_log_file: str, logger: Optional[logging.Logger] = None) -> bool:
    """
    Run KPCT with parallel processing and fallback to sequential if parallel fails.

    Parameters
    ----------
    kpct_input_file : str
        Path to the KPCT input file
    savedir : str
        Directory to save outputs
    kpct_outprefix : str
        Prefix for KPCT output files
    ncpus : int
        Number of CPU cores to use for parallel processing
    logger : Optional[logging.Logger], optional
        Logger instance for logging progress

    Returns
    -------
    bool
        True if KPCT ran successfully, False otherwise
    """

    if ncpus > 1:
        if logger:
            logger.info("Attempting parallel KPCT processing")

        parallel_success = run_kpct_parallel(kpct_input_file, savedir, kpct_outprefix, ncpus, logger)

        if parallel_success:
            return True
        else:
            if logger:
                logger.warning("Parallel KPCT processing failed, falling back to sequential processing")
            typer.secho("⚠️ [WARNING] Parallel processing failed, trying sequential approach", fg="yellow")


    if logger:
        logger.info("Running KPCT in sequential mode")

    return run_kpct(kpct_input_file, savedir, kpct_outprefix, resource_log_file, logger)


app = typer.Typer()

@app.command()
def pipeline(genomedir: str,
             savedir: str,
             ncpus: int=16,
             adapt_headers: bool=False,
             del_tmp: bool=True,
             calculate_complementarity: int=0,
             lowmem: bool = typer.Option(False, "--lowmem", help="Run emapper with reduced memory footprint, omitting --dbmem flag."),
             verbose: bool = typer.Option(False, "--verbose", help="Enable verbose output with detailed progress information."),
             ) -> None:
    """
    Run the ModuComp pipeline on a directory of genome files.

    This function executes the complete ModuComp workflow:
    1. Validates input genome files (.faa format)
    2. Optionally adapts FASTA headers or copies files to temporary directory
    3. Merges all genome files into a single file for annotation
    4. Runs eggNOG-mapper to annotate proteins with KO terms
    5. Creates a KO matrix from annotations
    6. Converts KO matrix to KPCT format
    7. Runs KPCT (KEGG Pathways Completeness Tool) with parallel processing
    8. Creates module completeness matrix
    9. Optionally generates complementarity reports for N-member combinations
    10. Cleans up temporary files if requested

    Parameters
    ----------
    genomedir : str
        Directory containing genome files in FAA format
    savedir : str
        Directory to save outputs
    ncpus : int, optional
        Number of CPUs to use for parallel processing, by default 16
    adapt_headers : bool, optional
        Whether to adapt FASTA headers to format: >genomeName|protein_N.
        Set to True if headers need standardization for eggNOG-mapper, by default False
    del_tmp : bool, optional
        Whether to delete temporary files after completion, by default True
    calculate_complementarity : int, optional
        Calculate complementarity between N genomes (0 = disable, 2+ = N-member combinations).
        When set to N (>=2), generates reports showing modules that are complete
        only in combinations of N genomes but not in individuals, by default 0
    lowmem : bool, optional
        Run emapper with reduced memory footprint by omitting --dbmem flag, by default False

    Raises
    ------
    SystemExit
        If required files are not found or processing steps fail
    """
    # Setup logging first to capture everything
    logger = setup_logging(savedir)

    # Setup resource monitoring
    resource_log_file = setup_resource_logging(savedir)
    if logger:
        logger.info(f"Resource monitoring enabled. Log file: {resource_log_file}")

    # Run the main pipeline logic
    _run_pipeline_core(genomedir, savedir, ncpus, adapt_headers, del_tmp,
                      calculate_complementarity, lowmem, verbose, logger, resource_log_file)


def _run_pipeline_core(genomedir: str, savedir: str, ncpus: int, adapt_headers: bool,
                      del_tmp: bool, calculate_complementarity: int, lowmem: bool,
                      verbose: bool, logger: logging.Logger, resource_log_file: str) -> None:
    """
    Core pipeline logic separated for resource monitoring.
    """
    start_time = time.time()

    greetings(verbose)
    conditional_output("\n📋 Initializing pipeline...", "green", verbose)

    # Convert to absolute paths
    genomedir = os.path.abspath(genomedir)
    savedir = os.path.abspath(savedir)

    # Create output directory
    create_output_dir(savedir, verbose)

    # Log pipeline parameters
    logger.info(f"Pipeline started with parameters:")
    logger.info(f"  - Genome directory: {genomedir}")
    logger.info(f"  - Output directory: {savedir}")
    logger.info(f"  - CPU cores: {ncpus}")
    logger.info(f"  - Adapt headers: {adapt_headers}")
    logger.info(f"  - Delete temporary files: {del_tmp}")
    logger.info(f"  - Calculate complementarity: {calculate_complementarity}")
    logger.info(f"  - Low memory mode: {lowmem}")

    # Check if all outputs already exist
    if check_final_reports_exist(savedir, calculate_complementarity, logger):
        conditional_output("✅ [OK] All output files already exist. Skipping processing.", "green", verbose)
        if not del_tmp:
            conditional_output("ℹ️ Keeping temporary files as requested.", "blue", verbose)
        logger.info("Pipeline skipped as all output files already exist")
        return

    # Validate input genomes
    conditional_output(f"🔍 Checking genome directory: {genomedir}", "yellow", verbose)
    how_many_genomes(genomedir, verbose)
    n_genomes = len(get_path_to_each_genome(genomedir))
    logger.info(f"Found {n_genomes} genome files in {genomedir}")

    # Create temporary directory
    create_tmp_dir(savedir, verbose)
    logger.info(f"Temporary directory created/verified")

    # Define file paths
    emapper_annotation_file = f"{savedir}/emapper_out.emapper.annotations"
    tmp_emapper_output_dir = f"{get_tmp_dir(savedir)}/emapper_output"
    tmp_emapper_file = f"{tmp_emapper_output_dir}/emapper_out.emapper.annotations"
    ko_matrix_path = f"{savedir}/kos_matrix.csv"

    # Process annotations and create KO matrix
    if os.path.exists(ko_matrix_path):
        logger.info(f"KO matrix already exists: {ko_matrix_path}")
        conditional_output(f"✅ [OK] Using existing KO matrix: {ko_matrix_path}", "white", verbose)
    else:
        # Check for existing emapper annotations
        if os.path.exists(emapper_annotation_file):
            logger.info(f"Emapper annotations already exist: {emapper_annotation_file}")
            conditional_output(f"✅ [OK] Using existing emapper annotations: {emapper_annotation_file}", "white", verbose)
        elif os.path.exists(tmp_emapper_file):
            logger.info(f"Emapper annotations found in temp directory: {tmp_emapper_file}")
            conditional_output(f"✅ [OK] Using existing emapper annotations from temp directory", "white", verbose)
            # Copy to final location
            try:
                shutil.copy(tmp_emapper_file, emapper_annotation_file)
                logger.info(f"Copied emapper annotations to: {emapper_annotation_file}")
            except Exception as e:
                logger.warning(f"Failed to copy emapper annotations: {str(e)}")
        else:
            # Need to run the full annotation pipeline
            # Prepare genome files
            if adapt_headers:
                logger.info("Starting to adapt fasta headers")
                adapt_fasta_headers(genomedir, savedir, verbose)
                logger.info("Completed adapting fasta headers")
            else:
                logger.info("Copying FAA files to temporary directory")
                copy_faa_to_tmp(genomedir, savedir, verbose)
                logger.info("Completed copying FAA files")

            # Merge genomes
            logger.info("Starting genome merging")
            merge_success = merge_genomes(savedir, logger, verbose)
            if not merge_success:
                logger.error("Failed to merge genomes. Exiting pipeline.")
                typer.secho("❌ [ERROR] Failed to merge genomes. Exiting pipeline.", fg="red")
                return

            # Run eggNOG-mapper
            logger.info(f"Starting eMapper with {ncpus} CPUs")
            emapper_success = run_emapper(savedir, ncpus, resource_log_file, lowmem, logger, verbose)
            if not emapper_success:
                logger.error("Failed to run emapper. Exiting pipeline.")
                typer.secho("❌ [ERROR] Failed to run emapper. Exiting pipeline.", fg="red")
                return

        # Create KO matrix from annotations
        logger.info(f"Creating KO matrix from eMapper annotations: {emapper_annotation_file}")
        create_ko_matrix_from_emapper_annotation(emapper_annotation_file, ko_matrix_path, logger)
        logger.info(f"KO matrix created: {ko_matrix_path}")

    # Process module completeness
    module_completeness_file = f"{savedir}/module_completeness.tsv"

    if os.path.exists(module_completeness_file):
        logger.info(f"Module completeness matrix already exists: {module_completeness_file}")
        typer.secho(f"✅ [OK] Using existing module completeness matrix: {module_completeness_file}", fg="white")
    else:
        # Set up KPCT processing
        kpct_outprefix = "output_give_completeness"
        kpct_input_file = os.path.join(savedir, "ko_file_for_kpct.txt")

        # Check if KPCT output already exists
        possible_kpct_files = [
            os.path.join(savedir, f"{kpct_outprefix}_contigs.with_weights.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_pathways.with_weights.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_contigs.tsv"),
            os.path.join(savedir, f"{kpct_outprefix}_pathways.tsv")
        ]
        kpct_file_exists = any(os.path.exists(f) for f in possible_kpct_files)

        # Convert KO matrix to KPCT format if needed
        if not os.path.exists(kpct_input_file):
            logger.info(f"Converting KO matrix to KPCT format: {ko_matrix_path}")
            ko_matrix_to_kpct_format(ko_matrix_path, savedir, calculate_complementarity, logger)
        else:
            logger.info(f"KPCT input file already exists: {kpct_input_file}")
            typer.secho(f"✅ [OK] Using existing KPCT input file: {kpct_input_file}", fg="white")

        # Run KPCT if needed
        if not kpct_file_exists:
            # Check KPCT installation
            if not check_kpct_installed(logger):
                return

            # Run KPCT with parallel processing
            logger.info(f"Running KPCT with parallel processing on file: {kpct_input_file}")
            kpct_success = run_kpct_with_fallback(kpct_input_file, savedir, kpct_outprefix, ncpus, resource_log_file, logger)
            if not kpct_success:
                return
        else:
            logger.info(f"KPCT output file(s) already exist with prefix '{kpct_outprefix}'")
            typer.secho(f"✅ [OK] Using existing KPCT output files with prefix '{kpct_outprefix}'", fg="white")

        # Create module completeness matrix
        logger.info(f"Creating module completeness matrix")
        create_module_completeness_matrix(savedir, kpct_outprefix, logger)

    # Generate complementarity reports if requested
    if calculate_complementarity >= 2:
        logger.info(f"Generating complementarity reports for up to {calculate_complementarity}-member combinations")

        # Generate reports for each combination size
        for n_members in range(2, calculate_complementarity + 1):
            complementarity_report_file = f"{savedir}/module_completeness_complementarity_{n_members}member.tsv"
            if os.path.exists(complementarity_report_file):
                logger.info(f"Complementarity report for {n_members}-member combinations already exists: {complementarity_report_file}")
                typer.secho(f"✅ [OK] Using existing {n_members}-member complementarity report: {complementarity_report_file}", fg="white")
            else:
                logger.info(f"Generating complementarity report for {n_members}-member combinations")
                generate_complementarity_report(savedir, n_members, logger)

    # Clean up temporary files if requested
    if del_tmp:
        logger.info("Cleaning up temporary files")
        remove_temp_files(savedir, logger)

    # Generate final resource usage summary
    log_final_resource_summary(resource_log_file, start_time, logger, verbose)

    end_time = time.time()
    elapsed_time = end_time - start_time
    typer.secho(f"\n✅ Pipeline completed successfully in {elapsed_time:.2f} seconds", fg="green", bold=True)
    logger.info(f"Pipeline completed in {elapsed_time:.2f} seconds")


@app.command()
def analyze_ko_matrix(
    kos_matrix: str,
    savedir: str,
    calculate_complementarity: int=0,
    kpct_outprefix: str="output_give_completeness",
    del_tmp: bool=True,
    ncpus: int=16,
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose output with detailed progress information."),
    ) -> None:
    """
    Run module completeness analysis on a pre-existing KO matrix file using the KEGG Pathways Completeness Tool (KPCT).

    This command provides an alternative entry point when you already have a KO matrix
    and want to skip the protein annotation steps. The workflow includes:
    1. Validates the input KO matrix file
    2. Converts KO matrix to KPCT input format
    3. Runs KPCT with parallel processing to calculate module completeness
    4. Creates module completeness matrix from KPCT output
    5. Optionally generates complementarity reports for N-member combinations
    6. Cleans up temporary files if requested

    The KO matrix should be in CSV/TSV format with:
    - First column: 'taxon_oid' containing genome identifiers
    - Remaining columns: KO identifiers (e.g., K00001, K00002, etc.)
    - Values: counts of proteins for each KO in each genome

    Parameters
    ----------
    kos_matrix : str
        Path to the KO matrix CSV/TSV file. Expected format: rows=genomes, columns=KO terms
    savedir : str
        Directory to save outputs and intermediate files
    calculate_complementarity : int, optional
        Calculate complementarity between N genomes (0 = disable, 2+ = N-member combinations).
        Generates reports showing modules complete only in combinations, by default 0
    kpct_outprefix : str, optional
        Prefix for KPCT output files, by default "output_give_completeness"
    del_tmp : bool, optional
        Whether to delete temporary files after completion, by default True
    ncpus : int, optional
        Number of CPUs to use for parallel KPCT processing, by default 16

    Raises
    ------
    SystemExit
        If KO matrix file is not found or KPCT processing fails
    """

    start_time = time.time()
    kos_matrix = os.path.abspath(kos_matrix)
    savedir = os.path.abspath(savedir)

    # Setup logging
    logger = setup_logging(savedir)

    # Setup resource monitoring
    resource_log_file = setup_resource_logging(savedir)
    if logger:
        logger.info(f"Resource monitoring enabled. Log file: {resource_log_file}")

    greetings(verbose)
    conditional_output("\n📋 Initializing KO matrix analysis...", "green", verbose)


    if not os.path.exists(kos_matrix):
        typer.secho(f"❌ [ERROR] KO matrix file not found at: {kos_matrix}", fg="red")
        exit(1)


    if not os.path.exists(savedir):
        os.makedirs(savedir)


    logger = setup_logging(savedir)

    if logger:
        logger.info("Starting moducomp KO matrix analysis")
        logger.info(f"KO matrix file: {kos_matrix}")
        logger.info(f"Output directory: {savedir}")
        logger.info(f"Calculate complementarity: {calculate_complementarity}")
        logger.info(f"KPCT output prefix: {kpct_outprefix}")
        logger.info(f"CPU cores: {ncpus}")


    if check_final_reports_exist(savedir, calculate_complementarity, logger):
        typer.secho("✅ [OK] All output files already exist. Skipping processing.", fg="green")
        logger.info("Analysis skipped as all output files already exist")
        return


    output_matrix_path = os.path.join(savedir, "kos_matrix.csv")
    if not os.path.exists(output_matrix_path):
        try:
            shutil.copy(kos_matrix, output_matrix_path)
            if logger:
                logger.info(f"Copied KO matrix to output directory: {output_matrix_path}")
        except Exception as e:
            if logger:
                logger.warning(f"Failed to copy KO matrix to output directory: {str(e)}")
    else:
        if logger:
            logger.info(f"KO matrix already exists in output directory: {output_matrix_path}")


    kpct_input_file = os.path.join(savedir, "ko_file_for_kpct.txt")
    module_completeness_file = os.path.join(savedir, "module_completeness.tsv")


    possible_kpct_files = [
        os.path.join(savedir, f"{kpct_outprefix}_contigs.with_weights.tsv"),
        os.path.join(savedir, f"{kpct_outprefix}_pathways.with_weights.tsv"),
        os.path.join(savedir, f"{kpct_outprefix}_contigs.tsv"),
        os.path.join(savedir, f"{kpct_outprefix}_pathways.tsv")
    ]
    kpct_file_exists = any(os.path.exists(f) for f in possible_kpct_files)

    try:

        if not os.path.exists(kpct_input_file):
            logger.info(f"Converting KO matrix to KPCT format: {kos_matrix}")
            ko_matrix_to_kpct_format(kos_matrix, savedir, calculate_complementarity, logger)
        else:
            logger.info(f"KPCT input file already exists: {kpct_input_file}")
            typer.secho(f"✅ [OK] Using existing KPCT input file: {kpct_input_file}", fg="white")


        if not kpct_file_exists:

            if not check_kpct_installed(logger):
                exit(1)


            logger.info(f"Running KPCT with parallel processing on file: {kpct_input_file}")
            kpct_success = run_kpct_with_fallback(kpct_input_file, savedir, kpct_outprefix, ncpus, resource_log_file, logger)
            if not kpct_success:
                exit(1)
        else:
            logger.info(f"KPCT output file(s) already exist with prefix '{kpct_outprefix}'")
            typer.secho(f"✅ [OK] Using existing KPCT output files with prefix '{kpct_outprefix}'", fg="white")


        if not os.path.exists(module_completeness_file):
            if logger:
                logger.info(f"Creating module completeness matrix")
            create_module_completeness_matrix(savedir, kpct_outprefix, logger)
        else:
            if logger:
                logger.info(f"Module completeness matrix already exists: {module_completeness_file}")
            typer.secho(f"✅ [OK] Using existing module completeness matrix: {module_completeness_file}", fg="white")


        if calculate_complementarity >= 2:
            if logger:
                logger.info(f"Generating complementarity reports for up to {calculate_complementarity}-member combinations")


            for n_members in range(2, calculate_complementarity + 1):
                complementarity_report_file = f"{savedir}/module_completeness_complementarity_{n_members}member.tsv"
                if os.path.exists(complementarity_report_file):
                    logger.info(f"Complementarity report for {n_members}-member combinations already exists: {complementarity_report_file}")
                    typer.secho(f"✅ [OK] Using existing {n_members}-member complementarity report: {complementarity_report_file}", fg="white")
                else:
                    logger.info(f"Generating complementarity report for {n_members}-member combinations")
                    generate_complementarity_report(savedir, n_members, logger)


        if del_tmp:
            if logger:
                logger.info("Cleaning up temporary files")
            remove_temp_files(savedir, logger)

        # Generate final resource usage summary
        log_final_resource_summary(resource_log_file, start_time, logger, verbose)

        end_time = time.time()

        typer.secho(f"\n✅ KO matrix KPCT analysis completed in {end_time - start_time:.2f} seconds", fg="green")
        typer.secho(f" Output directory: {savedir}", fg="green")

        if logger:
            logger.info(f"Analysis completed in {end_time - start_time:.2f} seconds")
            logger.info(f"Output directory: {savedir}")

    except Exception as e:
        if logger:
            logger.error(f"Error in KPCT analysis: {str(e)}", exc_info=True)
        typer.secho(f"❌ [ERROR] Error in KPCT analysis: {str(e)}", fg="red")
        exit(1)


if __name__ == "__main__":
    app()