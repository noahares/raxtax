import os
import argparse
import csv
import common

# Define main function to handle different sample sizes
def main(query, db, raxtax, sintax, thread_counts, repetitions, output_dir):

    csv_file = os.path.join(output_dir, "time_memory.csv")
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Threads', 'RuntimeSeconds', 'MaxMemoryBytes', 'Tool'])
    for num_threads in thread_counts:
        for i in range(repetitions):
            print(f"Running with {num_threads} threads, repetition {i+1}")
            # Prepare file paths
            raxtax_dir = os.path.join(output_dir, f"raxtax_{num_threads}_10pct_rep{i+1}")
            # sintax_dir = os.path.join(output_dir, f"sintax_{num_threads}_10pct_rep{i+1}")

            # Run the external program on the 90% file and measure performance
            r_runtime, r_max_memory = common.build_raxtax_command(raxtax, query, db, raxtax_dir, num_threads)
            s_runtime, s_max_memory = common.build_sintax_command(sintax, query, db, raxtax_dir, num_threads)
            print(f"RaxTax Runtime: {r_runtime:.2f} seconds, Max Memory: {r_max_memory / (1024 * 1024):.2f} MB")
            print(f"Sintax Runtime: {s_runtime:.2f} seconds, Max Memory: {s_max_memory / (1024 * 1024):.2f} MB")
            common.write_to_csv(csv_file, num_threads, r_runtime, r_max_memory, "RAxTax")
            common.write_to_csv(csv_file, num_threads, s_runtime, s_max_memory, "Sintax")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate speedup")
    parser.add_argument("-i", dest="input_fasta", type=str, help="Path to the input FASTA file")
    parser.add_argument("-d", dest="input_db", type=str, help="Path to the input FASTA db file")
    parser.add_argument("--raxtax", dest="raxtax", type=str, help="Path to raxtax")
    parser.add_argument("--sintax", dest="sintax", type=str, help="Path to usearch")
    parser.add_argument("-t", dest="thread_counts", nargs="+", type=int, default=[1, 2, 4, 8], help="List of thread counts")
    parser.add_argument("-r", dest="repetitions", type=int, default=3, help="Number of repetitions per sample size")
    parser.add_argument("-o", dest="output_dir", type=str, default="output_files", help="Directory to store output files")

    args = parser.parse_args()

    # Make sure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Run the main function
    main(args.input_fasta, args.input_db, args.raxtax, args.sintax, args.thread_counts, args.repetitions, args.output_dir)
