import os
import argparse
import csv
import common

# Define main function to handle different sample sizes
def main(input_fasta, raxtax, sintax, sample_sizes, repetitions, output_dir):

    csv_file = os.path.join(output_dir, "time_memory.csv")
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Rep', 'SampleSize', 'RuntimeSeconds', 'MaxMemoryBytes', 'Tool'])
    for num_samples in sample_sizes:
        for i in range(repetitions):
            print(f"Sampling {num_samples} sequences, repetition {i+1}")
            # Prepare file paths
            output_90 = os.path.join(output_dir, f"{num_samples}_90pct_rep{i+1}.fasta")
            output_10 = os.path.join(output_dir, f"{num_samples}_10pct_rep{i+1}.fasta")
            raxtax_dir = os.path.join(output_dir, f"raxtax_{num_samples}_10pct_rep{i+1}")

            # Sample sequences and split into 90% and 10%
            sample_90, sample_10 = common.sample_fasta(input_fasta, num_samples)
            with open(output_90, "w") as out90:
                SeqIO.write(sample_90, out90, "fasta")

            with open(output_10, "w") as out10:
                SeqIO.write(sample_10, out10, "fasta")

            # Run the external program on the 90% file and measure performance
            r_runtime, r_max_memory = common.build_raxtax_command(raxtax, output_10, output_90, raxtax_dir, 48)
            print(f"RaxTax Runtime: {r_runtime:.2f} seconds, Max Memory: {r_max_memory / (1024 * 1024):.2f} MB")
            common.write_to_csv(csv_file, i, num_samples, r_runtime, r_max_memory, "raxtax")
            s_runtime, s_max_memory = common.build_sintax_command(sintax, output_10, output_90, raxtax_dir, 48)
            print(f"Sintax Runtime: {s_runtime:.2f} seconds, Max Memory: {s_max_memory / (1024 * 1024):.2f} MB")
            common.write_to_csv(csv_file, i, num_samples, s_runtime, s_max_memory, "sintax")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample sequences from a FASTA file, split into 90%/10%, and run a program on them.")
    parser.add_argument("-i", dest="input_fasta", type=str, help="Path to the input FASTA file")
    parser.add_argument("--raxtax", dest="raxtax", type=str, help="Path to raxtax")
    parser.add_argument("--sintax", dest="sintax", type=str, help="Path to usearch")
    parser.add_argument("-s", dest="sample_sizes", nargs="+", type=int, default=[50000, 100000, 200000, 500000, 1000000], help="List of sample sizes")
    parser.add_argument("-r", dest="repetitions", type=int, default=3, help="Number of repetitions per sample size")
    parser.add_argument("-o", dest="output_dir", type=str, default="output_files", help="Directory to store output files")

    args = parser.parse_args()

    # Make sure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Run the main function
    main(args.input_fasta, args.raxtax, args.sintax, args.sample_sizes, args.repetitions, args.output_dir)
