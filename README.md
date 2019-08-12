A fast implementation of *de novo* clustering for error-prone full-length amplicon sequences

Includes two different approaches - hierarchical clustering of pairwise distances, and a locus-based "denoising" approach

Right now, only the hierarchical approach is implemented

Installation
------------

    git clone http://github.com/txje/flac
    cd flac
    mkdir incl
    cd incl
    git clone http://github.com/attractivechaos/klib
    cd ..
    make

To run a quick test:

    ./test.sh

Usage
-----

    ./flac -q <reads.fq> -r <ref.fa> -p <read_to_ref.paf>
    
    Usage: flac [options]
    Options:
      -q: FASTA/Q[.gz] file with reads
      -r: Reference FASTA/Q[.gz] or precomputed index file
      -p, --paf: PAF alignment file, or '-' for stdin
      -d, --distance-matrix: Write (or read, if it exists) distance matrix file
      -b, --begin: begin position to evaluate (in single ref seq)
      -e, --end: end position to evaluate (in single ref seq)

      Alignment parameters:
        -t: Threads (default: 1)
        -s: Sequence type preset
            map-ont: Oxford Nanopore (default)
            map-pb:  Pacbio
        -f, --align-fraction: Portion of a read that must align properly (defaults to --align-length threshold)
        -l, --align-length: Minimum aligned bp (default: 100)
        -a, --align-accuracy: Minimum accuracy of aligned portion of a read (default: 0.6)

      -v, --verbose: verbose
      -h, --help: show this
      --version: show version information
