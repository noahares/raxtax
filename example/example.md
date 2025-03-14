# Diptera COX1 Example

The files `diptera_references.fasta` and `diptera_queries.fasta` contain *Diptera* sequences for a quick example run of `raxtax`.

From the project root (otherwise adjust the paths) run:
```sh
cargo run --profile=ultra -- -d example/diptera_references.fasta -i example/diptera_queries.fasta -o example/example_run
```

This creates a new folder `example/example_run` with the taxonomic assignments and confidence values for each query in `raxtax.out` and various log messages (including exact sequence matches) in `raxtax.log`.
