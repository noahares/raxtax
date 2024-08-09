# RAxTax - RAxTax Accelerates Taxonomic Classification

This project is heavily inspired by the SINTAX algorithm [[1]](#1).

## Usage

For maximum performance, build the program with `cargo build --profile=ultra`.

```sh
Usage: raxtax [OPTIONS] --database-path <DATABASE_PATH> --query-file <QUERY_FILE>

Options:
  -d, --database-path <DATABASE_PATH>  Path to the database fasta or bin file
  -i, --query-file <QUERY_FILE>        Path to the query file
      --skip-exact-matches             If used for mislabling analysis, you want to skip exact sequence matches
      --tsv                            Output primary result file in tsv format
      --make-db                        Create a binary database to load instead of a fasta file for repeated execution
  -t, --threads <THREADS>              Number of threads
                                       If 0, uses all available threads [default: 0]
  -o, --prefix <PREFIX>                Output prefix
      --redo                           Force override of existing output files
  -v, --verbose...                     Increase logging verbosity
  -q, --quiet...                       Decrease logging verbosity
  -h, --help                           Print help
  -V, --version                        Print version
```

## Format

Input files may be provided as gzip compressed archives.

### Input Database

The input format for the database file is FASTA. 
Sequence identifier should have the form `tax=<lineage>;`. 
Everything after `tax=` is parsed as a comma-separated list of lineage nodes and is terminated by a semicolon.
Lineages may have different depth, the only requirement is that they can be parsed into a multi-furcating tree.
We use phylum to sequence for the examples in this README to aid readability.
For example, an entry may look like this:

```sh
# example sequence
>metadata;tax=phylum,class,order,family,genus,species;
ACTCGATAC
```

### Input Query

The format for query sequences is also FASTA, but more relaxed than the database format:

```sh
# example sequence
>label
ACTCGATAC
```

### Output

RAxTax will produce 2 primary output files under the prefix specified with `-o` (defaults to `name_of_query_file.out/`).

1. `<PREFIX>/raxtax.out` is the full result of the analysis. It contains for each query sequence a line for each database sequence where the confidence value is above 0.01 (confidence values are between 0 and 1).
**If no database sequence fullfills this criterion, a single line containing the best match is printed**.
In this case, values are rounded up to 0.01.
The format is (tab separated):

```sh
# regular output
<query_label>    <p>,<c>,<o>,<f>,<g>,<s>    <c_p>,<c_c>,<c_o>,<c_f>,<c_g>,<c_s>  <local-signal>  <global-signal>
```

The first part is simply the query label.
The second part is the taxonomical lineage of the respective database sequence.
The third part contains the confidence values for each level of the taxonomical lineage.
**It is important to understand that these values are always relative to the sequences in the database and therefore should be interpreted carefully.**
To this end, we include a fourth and fifth value indicating the confidence in the reported lineage (local signal) and confidence in the confidence values themselves on sequence level (global signal).
These are again between 0 and 1, where 1 indicates high confidence.
For more information, see the manuscript.

2. `<PREFIX>/raxtax.log` is the log file where more or less useful information accumulates.
With the default command line parameters, only warnings and errors will be collected.
With `-v` additional information about runtime and the size of the database are printed.
With `-vv` debug messages are also included.
Generally, if a warning or error occurs, the program will inform you through `stderr` and refer you to the log file if needed.
This file also contains information about exact matches and inconsistent lineages (possible mislabeling).

4. (_optional_ via `--tsv`) `<PREFIX>/raxtax.tsv` is pretty much the same as the first output file but slightly more convenient for viewing in your favorite spreadsheet editor.
In this file, the taxonomical lineage and confidence values are interleaved and the query sequence is also printed at the end:

```sh
<query_label>    <p>    <c_p>   <c> <c_c>   <o> <c_o>   <f> <c_f>   <g> <c_g>   <s> <c_s>  <local-signal>  <global-signal>   <query_sequence>
```

## Other Options

`--skip-exact-matches` may be useful when running the database against itself to identify mislabeled sequences. Per default, RAxTax skips over exact sequences matches if there is **exactly one match** and outputs a confidence of 1.0 for the exact match.
This option makes it so that any exact match is not considered for the analysis of a query sequence.

`--make-db` can be used if you want to run the program with the same reference database for many different query files.
If the reference database is large this will save significant time on repeat execution.

`--threads` may be ommited most of the time and RAxTax will use as many cores as your system has available. Because the analysis is _embarrassingly parallel_, this is a sensible default.
However, if you experience problems due to hyper-threading, you might want to reduce the number of threads, to increase parallel efficiency.

`--redo` will enable overwriting of existing output files. **Use at your own risk!**

## Important Implementation Details

We suggest a threshold of 0.01 for confidence values to be considered (`F64_OUTPUT_ACCURACY` also in `src/utils.rs`).
For technical reasons this is the number of digits after the decimal point, so currently this is 2.

If the database contains duplicate sequences that have different lineages above the lowest taxonomic level a warning will be emitted.

## Why does this exist if there are so many taxonomical classifiers already, and how does it work?

We will soon publish a manuscript about this method and what we use it for.

## References
<a id="1">[1]</a>
Edgar, Robert C. "SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences." biorxiv (2016): 074161.

## Copyright
 This work is licensed under CC BY-NC-SA 4.0. 
 To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
