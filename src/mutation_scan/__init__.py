"""MutationScan: protein variant calling and descriptive mutation co-occurrence.

Stages
------
1. ``references``   discover target genes from the reference proteins supplied
2. ``extract``      tblastn translating alignment: genome DNA -> target proteins
3. ``variants``     global alignment + gap-aware residue counter -> substitutions
4. ``cooccurrence`` per-genome counts of which genes are mutated together
5. ``epistasis``    genome mutation maps + statistical epistasis networks

No gene, genome, organism, drug or structure is hardcoded anywhere in this
package. The targets of a run are whatever reference FASTA files the user
supplies, and every output is a count or a measurement -- the tool reports what
is in the data and interprets none of it.
"""

# Must happen before any submodule imports numpy or scipy: BLAS reads these
# once, at load time, so pinning later is silently a no-op. This is why the
# call sits here in the package root rather than in ``mutation_scan.parallel``.
# Values the user set explicitly are left alone -- this only fills in blanks.
from .parallel import pin_blas_threads as _pin_blas_threads

_pin_blas_threads()

__version__ = "1.0.0"

__all__ = ["__version__"]
