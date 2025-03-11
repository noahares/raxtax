import os
import argparse
import csv
import common
from Bio import SeqIO

# Define main function to handle different sample sizes
def main(input_fasta, raxtax, thread_counts, num_samples, repetitions, weak_scaling, fixed_query, output_dir):

    csv_file = os.path.join(output_dir, "time_memory.csv")
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Rep', 'Threads', 'RuntimeSeconds', 'MaxMemoryBytes', 'Tool'])
    for i in range(repetitions):
        output_90 = os.path.join(output_dir, f"90pct_rep{i+1}.fasta")
        factr = max(thread_counts) if weak_scaling else 1
        sample_90, sample_10 = common.sample_fasta(input_fasta, num_samples, factr, fixed_query)
        with open(output_90, "w") as out90:
            SeqIO.write(sample_90, out90, "fasta")
        if weak_scaling:
            for num_threads in thread_counts:
                output_10 = os.path.join(output_dir, f"{num_threads}_10pct_rep{i+1}.fasta")
                with open(output_10, "w") as out10:
                    SeqIO.write(sample_10[:((num_threads) * 2000)], out10, "fasta")
        else:
            output_10 = os.path.join(output_dir, f"query_rep{i+1}.fasta")
            with open(output_10, "w") as out10:
                SeqIO.write(sample_10, out10, "fasta")
    for num_threads in thread_counts:
        for i in range(repetitions):
            print(f"Running with {num_threads} threads, repetition {i+1}")
            # Prepare file paths
            output_10 = os.path.join(output_dir, f"{num_threads}_10pct_rep{i+1}.fasta") if weak_scaling else os.path.join(output_dir, f"query_rep{i+1}.fasta")
            output_90 = os.path.join(output_dir, f"90pct_rep{i+1}.fasta")
            raxtax_dir = os.path.join(output_dir, f"raxtax_{num_threads}_10pct_rep{i+1}")


            # Run the external program on the 90% file and measure performance
            r_runtime, r_max_memory = common.build_raxtax_command(raxtax, output_10, output_90, raxtax_dir, num_threads)
            print(f"RaxTax Runtime: {r_runtime:.2f} seconds, Max Memory: {r_max_memory / (1024 * 1024):.2f} MB")
            common.write_to_csv(csv_file, i, num_threads, r_runtime, r_max_memory, "raxtax")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate speedup")
    parser.add_argument("-i", dest="input_fasta", type=str, help="Path to the input FASTA file")
    parser.add_argument("--raxtax", dest="raxtax", type=str, help="Path to raxtax")
    parser.add_argument("-t", dest="thread_counts", nargs="+", type=int, default=[1, 2, 4, 8], help="List of thread counts")
    parser.add_argument("-s", dest="sample_size", type=int, default=100000, help="Sample size")
    parser.add_argument("-r", dest="repetitions", type=int, default=3, help="Number of repetitions per sample size")
    parser.add_argument("-o", dest="output_dir", type=str, default="output_files", help="Directory to store output files")
    parser.add_argument("-f", dest="fixed_query", type=bool, default=False, help="If true, keep query size fixed at 2000")
    parser.add_argument("-w", dest="weak_scaling", type=bool, default=False, help="If true, each thread adds 2000 queries")

    args = parser.parse_args()

    # Make sure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Run the main function
    main(args.input_fasta, args.raxtax, args.thread_counts, args.sample_size, args.repetitions, args.weak_scaling, args.fixed_query, args.output_dir)
