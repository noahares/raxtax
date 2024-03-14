# make-sintax-great-again

This is an implementation of the SINTAX algorithm [[1]](#1).

## Usage

```sh
Usage: make-sintax-great-again [OPTIONS] --sequence-file <SEQUENCE_FILE> --query-file <QUERY_FILE>

Options:
  -s, --sequence-file <SEQUENCE_FILE>      Path to the sequence file
  -d, --database-path <DATABASE_PATH>      Path to the existing database
  -z, --database-output <DATABASE_OUTPUT>  Path to the database output
  -i, --query-file <QUERY_FILE>            Path to the query file
  -n, --num-rounds <NUM_ROUNDS>            Number of rounds per query [default: 100]
  -k, --num-k-mers <NUM_K_MERS>            Number of 8-mers [default: 32]
  -p, --threshold <THRESHOLD>              8-mer hit-threshold [default: 0.3333333333333333]
  -r, --num-results <NUM_RESULTS>          Number of output species per query [default: 5]
  -t, --threads <THREADS>                  Number of threads [default: 0]
      --seed <SEED>                        Seed [default: 42]
  -o, --output <OUTPUT>                    Output path
  -v, --verbose...                         Increase logging verbosity
  -q, --quiet...                           Decrease logging verbosity
  -h, --help                               Print help
  -V, --version                            Print version
```

## References
<a id="1">[1]</a>
Edgar, Robert C. "SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences." biorxiv (2016): 074161.
