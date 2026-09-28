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

__version__ = "1.0.0"

__all__ = ["__version__"]
