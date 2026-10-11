"""Benchmarking harness for MutationScan.

Two kinds of benchmark live here:

* :mod:`parallel_bench` -- the executor scaling
  benchmark (standalone script, run directly).
* The AMR tool comparison -- :mod:`run_benchmark`
  scores MutationScan against AMRFinderPlus,
  Mykrobe, Kleborate and StarAMR on a shared
  ground truth. Run it as a module so the
  package-relative imports resolve::

      python -m benchmarks.run_benchmark --help
"""
