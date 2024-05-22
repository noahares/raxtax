# RAxTax - RAxTax Accelerates Taxonomic Classification

This project is heavily inspired by the SINTAX algorithm [[1]](#1).

## Usage

For maximum performance, build the program with `cargo build --profile=ultra`.

```sh
Usage: raxtax [OPTIONS] --database-path <DATABASE_PATH> --query-file <QUERY_FILE>

Options:
  -d, --database-path <DATABASE_PATH>
          Path to the database fasta file
  -i, --query-file <QUERY_FILE>
          Path to the query file
  -c, --min-confidence <MIN_CONFIDENCE>
          Confidence threshold [default: 0.8]
      --skip-exact-matches
          If used for mislabling analysis, you want to skip exact sequence matches
      --tsv
          Output primary result file in tsv format
  -t, --threads <THREADS>
          Number of threads
          If 0, uses all available threads [default: 0]
  -o, --prefix <PREFIX>
          Output prefix
      --redo
          Force override of existing output files
  -v, --verbose...
          Increase logging verbosity
  -q, --quiet...
          Decrease logging verbosity
  -h, --help
          Print help
  -V, --version
          Print version
```

## Format

Input files may be provided as gzip compressed archives.

### Input Database

The input format for the database file is FASTA and follows the same convention as the input for SINTAX:

```sh
# example sequence
>metadata;tax=p:phylum,c:class,o:order,f:family,g:genus,s:species;
ACTCGATAC
```

Taxonomical lineages must be consistent. These labels lead to an error because `genus` belongs to both, `family` and `family2`.

```sh
>metadata;tax=p:phylum,c:class,o:order,f:family,g:genus,s:species;
ACTCGATAC
>metadata;tax=p:phylum,c:class,o:order,f:family2,g:genus,s:species2;
AATCGATTC
```

Taxonomical identifiers may not contain any characters out of the set \[,:;|\] after `tax=`.

### Input Query

The format for query sequences is also FASTA, but more relaxed than the database format:

```sh
# example sequence
>label
ACTCGATAC
```

### Output

RAxTax will produce 3 primary output files under the prefix specified with `-o` (defaults to `name_of_query_file.out/`).

1. `<PREFIX>/raxtax.out` is the full result of the analysis. It contains for each query sequence a line for each database sequence where the confidence value is above 0.01 (confidence values are between 0 and 1).
**If no database sequence fullfills this criterion, a single line containing `NA` will be printed.
The format is (tab separated):

```sh
# regular output
<query_label>    <p>|<c>|<o>|<f>|<g>|<s>    <c_p>|<c_c>|<c_o>|<c_f>|<c_g>|<c_s>
# no match above threshold
<query_label>   NA  NA
```

The first part is simply the query label.
The second part is the taxonomical lineage of the respective database sequence separated by '|'.
The last part contains the confidence values for each level of the taxonomical lineage.
**It is important to understand that these values are always relative to the sequences in the database and therefore should be interpreted carefully.**

2. `<PREFIX>/raxtax.confidence<c>.out` is a slightly truncated version of the first file where the database lineage and confidence values are discontinued if the confidence is bellow the threshold provided by `-c` (defaults to 0.8).

3. `<PREFIX>/raxtax.log` is the log file where more or less useful information accumulates.
With the default command line parameters, only warnings and errors will be collected.
With `-v` additional information about runtime and the size of the database are printed.
With `-vv` debug messages are also included.
Generally, if a warning or error occurs, the program will inform you through `stderr` and refer you to the log file if needed.

4. (_optional_ via `--tsv`) `<PREFIX>/raxtax.tsv` is pretty much the same as the first output file but slightly more convenient for viewing in your favorite spreadsheet editor.
In this file, the taxonomical lineage and confidence values are interleaved and the query sequence is also printed at the end:

```sh
<query_label>    <p>    <c_p>   <c> <c_c>   <o> <c_o>   <f> <c_f>   <g> <c_g>   <s> <c_s>   <query_sequence>
```

## Other Options

`--skip-exact-matches` may be useful when running the database against itself to identify mislabeled sequences. Per default, RAxTax skips over exact sequences matches and outputs a confidence of 1.0 for the exact match (or shared among multiple exact matches).
This option makes it so that any exact match is not considered for the analysis of a query sequence.

`--threads` may be ommited most of the time and RAxTax will use as many cores as your system has available. Because the analysis is _embarrassingly parallel_, this is a sensible default.
However, if you experience problems due to hyper-threading, you might want to reduce the number of threads, to increase parallel efficiency.

`--redo` will enable overwriting of existing output files. **Use at your own risk!**

## Important Implementation Details

We suggest that any class confidence value below 0.8 is meaningless, so results not satisfying this condition will be discarded.
If you disagree, you can change the value of `MIN_CLASS_CONFIDENCE` in `src/utils.rs` (and tell us your use-case because we are curious).

The same goes for the threshold of 0.01 for confidence values to be considered (`F64_OUTPUT_ACCURACY` also in `src/utils.rs`).
For technical reasons this is the number of digits after the decimal point, so currently this is 2.

If the database contains duplicate sequences that have different lineages above the species level, the confidence values will be wrong. A warning will be emitted.

## Why does this exist if there are so many taxonomical classifiers already, and how does it work?

We will soon publish a manuscript about this method and what we use it for.

## References
<a id="1">[1]</a>
Edgar, Robert C. "SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences." biorxiv (2016): 074161.
