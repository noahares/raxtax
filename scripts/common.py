import random
from Bio import SeqIO
import os
import subprocess
import time
import psutil
import csv

# Define a function to randomly sample sequences from a FASTA file
def sample_fasta(input_fasta, num_samples, output_90, output_10, fixed_query=False):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    sampled_sequences = random.sample(sequences, num_samples)

    random.shuffle(sampled_sequences)
    # Split the sampled sequences into 90% and 10%
    if fixed_query:
        split_point = int(len(sampled_sequences) - 2000)
    else:
        split_point = int(0.9 * len(sampled_sequences))
    sample_90 = sampled_sequences[:split_point]
    sample_10 = sampled_sequences[split_point:]

    # Write the 90% set to a file
    with open(output_90, "w") as out90:
        SeqIO.write(sample_90, out90, "fasta")

    # Write the 10% set to a file
    with open(output_10, "w") as out10:
        SeqIO.write(sample_10, out10, "fasta")

def build_raxtax_command(program, input_file, database_file, dir, num_threads):
    command = [program, "-i", input_file, "-d", database_file, "-o", dir, "-t", str(num_threads)]
    return run_program(command)

def build_sintax_command(program, input_file, database_file, dir, num_threads):
    db_path = os.path.splitext(database_file)[0] + ".udb"
    out_path = os.path.join(dir, "sintax.out")
    command_1 = [program, "--makeudb_usearch", database_file, "--output", db_path]
    command_2 = [program, "--sintax", input_file, "--db", db_path, "--tabbedout", out_path, "--threads", str(num_threads)]
    db_time, db_mem = run_program(command_1)
    exe_time, exe_mem = run_program(command_2)
    return db_time + exe_time, max(db_mem, exe_mem)

# Define a function to run the external program and measure runtime and memory usage
def run_program(command):
    # Start timing and memory usage monitoring
    start_time = time.time()
    process = psutil.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Monitor memory usage
    memory_usage = []
    try:
        while process.poll() is None:
            mem_info = process.memory_info().rss
            memory_usage.append(mem_info)
            time.sleep(0.1)  # Check every 0.1 seconds
    except psutil.NoSuchProcess:
        pass

    # Wait for process to finish and capture runtime
    stdout, stderr = process.communicate()
    end_time = time.time()

    runtime = end_time - start_time
    max_memory = max(memory_usage) if memory_usage else 0

    return runtime, max_memory
def write_to_csv(csv_file, field, runtime, memory, program):
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([field, runtime, memory, program])
