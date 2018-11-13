#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "incl/klib/kseq.h"
#include "hierarch.h"
#include "paf.h"

KSEQ_INIT(gzFile, gzread)

void version() {
  printf("FLAC (Full-Length Amplicon Clustering) version 0.1\n");
}

void usage() {
  printf("Usage: flac [options]\n");
  printf("Options:\n");
  printf("  -q: FASTA/Q[.gz] file with reads\n");
  printf("  -r: Reference FASTA/Q[.gz] or precomputed index file\n");
  printf("  --paf: PAF alignment file, or '-' for stdin\n");
  printf("  -t: Threads (default: 1)\n");
  printf("  -p: Sequence type preset\n");
  printf("      map-ont: Oxford Nanopore (default)\n");
  printf("      map-pb:  Pacbio\n");
  printf("  -f, --align-fraction: Portion of a read that must align properly (defaults to --align-length threshold)\n");
  printf("  -l, --align-length: Minimum aligned bp (default: 100)\n");
  printf("  -a, --align-accuracy: Minimum accuracy of aligned portion of a read (default: 0.6)\n");
  printf("  -v, --verbose: verbose\n");
  printf("  -h, --help: show this\n");
  printf("  --version: show version information\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "align-fraction",         required_argument, 0, 'f' },
  { "align-length",           required_argument, 0, 'l' },
  { "align-accuracy",         required_argument, 0, 'a' },
  { "verbose",                no_argument,       0, 'v' },
  { "help",                   no_argument,       0, 'h' },
  { "paf",                    required_argument, 0, 0 },
  { "version",                no_argument,       0, 0 },
  { 0, 0, 0, 0}
};

int parse_paf(char* paf_file) {
  fprintf(stderr, "Reading from alignment file '%s'\n", paf_file);
  paf_file_t *p = paf_open(paf_file);
  paf_rec_t r;
  int ret = paf_read(p, &r);
  while(ret == 0) {
    fprintf(af, "%s\t%d\t%d\t%d\t%c\t", r.qn, r.ql, r.qs, r.qe, "+-"[r.rev]);
    fprintf(af, "%s\t%d\t%d\t%d\t%d\t%d\t%d", r.tn, r.tl, r.ts, r.te, r.ml, r.bl, r.mq);
    fprintf(af, "\n");
    ret = paf_read(p, &r);
  }
  paf_close(p);
  fprintf(stderr, "Closing file '%s'\n", paf_file);
}

int main(int argc, char *argv[]) {
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;

  char* read_fasta = NULL;
  char* ref_fasta = NULL;
  char* paf_file = NULL;
  char* preset = "map-ont";
  float align_fraction = -1;
  float align_accuracy = 0.6;
  int align_length = 100;
  int n_threads = 1;
  int verbose = 0;

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "q:r:t:p:f:l:a:vh", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'q':
        read_fasta = optarg;
        break;
      case 'r':
        ref_fasta = optarg;
        break;
      case 't':
        n_threads = atoi(optarg);
        break;
      case 'p':
        preset = optarg;
        break;
      case 'f':
        align_fraction = atof(optarg);
        break;
      case 'a':
        align_accuracy = atof(optarg);
        break;
      case 'l':
        align_length = atoi(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'c':
        careful = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'q' || optopt == 'r' || optopt == 't' || optopt == 'f' || optopt == 'a' || optopt == 'p' || optopt == 'l')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        // as long as all the long arguments have characters too, I don't think this section will be used
        if (long_idx == 0) align_fraction = atof(optarg); // --align-fraction
        else if (long_idx == 1) align_length = 100; // --align-length
        else if (long_idx == 2) align_accuracy = atof(optarg); // --align-accuracy
        else if (long_idx == 3) verbose = 1; // --verbose
        else if (long_idx == 4) {usage(); return 0;} // --help
        else if (long_idx == 5) paf_file = optarg; // --paf
        else if (long_idx == 6) {version(); return 0;} // --version
        break;
      default:
        usage();
        return 1;
    }
  }

  // tests
  uint32_t n_items = 5;
  float **dists = malloc(sizeof(float*) * n_items);
  int i;
  for(i = 0; i < n_items; i++)
    dists[i] = malloc(sizeof(float) * n_items);
  // test case for points on a line like 01--2-34
  dists[0][0] = 0;
  dists[0][1] = 1;
  dists[0][2] = 4;
  dists[0][3] = 6;
  dists[0][4] = 7;

  dists[1][0] = 1;
  dists[1][1] = 0;
  dists[1][2] = 3;
  dists[1][3] = 5;
  dists[1][4] = 6;

  dists[2][0] = 4;
  dists[2][1] = 3;
  dists[2][2] = 0;
  dists[2][3] = 2;
  dists[2][4] = 3;

  dists[3][0] = 6;
  dists[3][1] = 5;
  dists[3][2] = 2;
  dists[3][3] = 0;
  dists[3][4] = 1;

  dists[4][0] = 7;
  dists[4][1] = 6;
  dists[4][2] = 3;
  dists[4][3] = 1;
  dists[4][4] = 0;
  agglomerate(dists, n_items);
}
