# sintax-ng

This is a modified implementation of the SINTAX algorithm [[1]](#1).

## Usage

```sh
Usage: sintax-ng [OPTIONS] --database-path <DATABASE_PATH> --query-file <QUERY_FILE>

Options:
  -d, --database-path <DATABASE_PATH>
          Path to the database fasta file
  -i, --query-file <QUERY_FILE>
          Path to the query file
  -k, --num-k-mers <NUM_K_MERS>
          Number of 8-mers [default: 32]
  -f, --min-hit-fraction <MIN_HIT_FRACTION>
          8-mer hit-threshold [default: 0.3333333333333333]
  -c, --min-confidence <MIN_CONFIDENCE>
          Confidence threshold [default: 0.8]
  -m, --max-target-seqs <MAX_TARGET_SEQS>
          Number of output species per query [default: 5]
      --skip-exact-matches
          If used for mislabling analysis, you want to skip exact sequence matches
  -t, --threads <THREADS>
          Number of threads
          If 0, uses all available threads [default: 0]
  -s, --seed <SEED>
          Seed [default: 42]
  -o, --output <OUTPUT>
          Output path
  -u, --confidence-output <CONFIDENCE_OUTPUT>
          confidence output path
  -l, --log-output <LOG_OUTPUT>
          log output path
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
