# RAxTax - RAxTax Accelerates Taxonomic Classification

This project is heavily inspired by the SINTAX algorithm [[1]](#1).

## Usage

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
  -s, --seed <SEED>
          Seed [default: 42]
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

## References
<a id="1">[1]</a>
Edgar, Robert C. "SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences." biorxiv (2016): 074161.
