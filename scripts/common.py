import os
import subprocess
import time
import psutil
import csv

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
