#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "incl/minimap2/minimap.h"
#include "paf.h"
#include "incl/klib/khash.h"
#include "incl/klib/kvec.h"
#include "hierarch.h"

#ifndef _kseq_
#define _kseq_

#include "incl/klib/kseq.h"

// init kseq struct
KSEQ_INIT(gzFile, gzread)

#endif

typedef struct {
  char* s;
  uint32_t l;
} fa_seq;

// map seq name to index into *seq vec
KHASH_MAP_INIT_STR(faHash, uint32_t)

void version() {
  printf("FLAC (Full-Length Amplicon Clustering) version 0.1\n");
}

void usage() {
  printf("Usage: flac [options]\n");
  printf("Options:\n");
  printf("  -q: FASTA/Q[.gz] file with reads\n");
  printf("  -r: Reference FASTA/Q[.gz] or precomputed index file\n");
  printf("  -p, --paf: PAF alignment file, or '-' for stdin\n");
  printf("  -t: Threads (default: 1)\n");
  printf("  -s: Sequence type preset\n");
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
  { "paf",                    required_argument, 0, 'p' },
  { "version",                no_argument,       0, 0 },
  { 0, 0, 0, 0}
};

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

  // ---------- options ----------
  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "q:r:t:p:s:f:l:a:vh", long_options, &long_idx)) != -1) {
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
      case 's':
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
      case 'p':
        paf_file = optarg;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'q' || optopt == 'r' || optopt == 't' || optopt == 'f' || optopt == 'a' || optopt == 'p' || optopt == 'l' || optopt == 's')
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

  if(read_fasta == NULL) {
    fprintf(stderr, "-q read FASTA is required\n");
    usage();
    return 1;
  }

  if(ref_fasta == NULL) {
    fprintf(stderr, "-r ref FASTA is required\n");
    usage();
    return 1;
  }

  // ---------- general initialization ----------
  int i, j, l;

  // capitalize standard nucleotides, turn anything else into N
  char norm[256];
  for(i = 0; i < 256; i++) {
    norm[i] = 'N';
  }
  norm[65] = 'A';
  norm[67] = 'C';
  norm[71] = 'G';
  norm[84] = 'T';
  norm[97] = 'A';
  norm[99] = 'C';
  norm[103] = 'G';
  norm[116] = 'T';

  char compl[256];
  for(i = 0; i < 256; i++) {
    compl[i] = 'N';
  }
  compl[65] = 'T';
  compl[67] = 'G';
  compl[71] = 'C';
  compl[84] = 'A';
  compl[97] = 'T';
  compl[99] = 'G';
  compl[103] = 'C';
  compl[116] = 'A';

  // ---------- load ref FASTA file ----------
  gzFile f = gzopen(ref_fasta, "r");
  kseq_t* seq = kseq_init(f);
  fprintf(stderr, "Loading ref FASTA file: %s\n", ref_fasta);

  khash_t(faHash) *refmap = kh_init(faHash);
  kvec_t(fa_seq) refs;
  kv_init(refs);

  khint_t bin; // hash bin (result of kh_put)
  int absent;
  char* a;
  fa_seq fs;
  fa_seq fn;

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    // make <seq>
    fs.s = malloc(sizeof(char)*l);
    memcpy(fs.s, seq->seq.s, sizeof(char)*l);
    fs.l = l;

    // add <seq> to refs vector
    kv_push(fa_seq, refs, fs);

    // add <name>:<refs idx> to refmap
    a = malloc(strlen(seq->name.s)+1);
    strcpy(a, seq->name.s);
    bin = kh_put(faHash, refmap, a, &absent);
    kh_val(refmap, bin) = kv_size(refs)-1;
  }

  fprintf(stderr, "Loaded %d sequences.\n", kv_size(refs));

  kseq_destroy(seq);
  gzclose(f);

  // ---------- load read FASTA file ----------
  f = gzopen(read_fasta, "r");
  seq = kseq_init(f);
  fprintf(stderr, "Loading read FASTA file: %s\n", read_fasta);

  khash_t(faHash) *readmap = kh_init(faHash);
  kvec_t(fa_seq) reads;
  kv_init(reads);
  kvec_t(fa_seq) readnames;
  kv_init(readnames);

  while ((l = kseq_read(seq)) >= 0) {
    // name: seq->name.s, seq: seq->seq.s, length: l
    //printf("Reading %s (%i bp).\n", seq->name.s, l);

    // make <seq>
    fs.s = malloc(sizeof(char)*l);
    memcpy(fs.s, seq->seq.s, sizeof(char)*l);
    fs.l = l;

    // add <seq> to reads vector
    kv_push(fa_seq, reads, fs);

    // add <name>:<reads idx> to readmap
    fn.s = malloc(strlen(seq->name.s)+1);
    strcpy(fn.s, seq->name.s);
    fn.l = strlen(fn.s);
    kv_push(fa_seq, readnames, fn);
    bin = kh_put(faHash, readmap, fn.s, &absent);
    kh_val(readmap, bin) = kv_size(reads)-1;
  }

  fprintf(stderr, "Loaded %d sequences.\n", kv_size(reads));

  kseq_destroy(seq);
  gzclose(f);

  // ---------- process from PAF file ----------
  int st = 2250;
  int en = 4450;
  // for simulated data:
  st=1100;
  en=2900;

  /*
   * matrix of reads x pileup, where uint8_t values are characters in this set:
   * [space]: not aligned (ascii 32)
   * -: deletion (ascii 45)
   * A, C, G, T (ascii 65, 67, 71, 84)
   * N (ascii 78, any other character in the input read is also turned into N)
   */
  uint8_t **matrix = malloc(sizeof(uint8_t*) * kv_size(reads));
  for(i = 0; i < kv_size(reads); i++) {
    matrix[i] = malloc((en-st) * sizeof(uint8_t));
    for(j = 0; j < (en-st); j++) matrix[i][j] = 32; // [space]
  }

  if(paf_file != NULL) {
    fprintf(stderr, "Reading from alignment file '%s'\n", paf_file);
    paf_file_t *p = paf_open(paf_file);
    paf_rec_t r;
    int ret = paf_read(p, &r);
    uint32_t val; // cigar op value
    uint32_t q; // position in query (read)
    uint32_t t; // position in target (ref)
    uint32_t qi; // query index in read list
    uint32_t ti; // target index into ref list
    char* qseq; // will point to the ENTIRE read seq
    char* tseq; // will point to the ENTIRE ref seq
    while(ret == 0) {
      /*
      fprintf(stdout, "%s\t%d\t%d\t%d\t%c\t", r.qn, r.ql, r.qs, r.qe, "+-"[r.rev]);
      fprintf(stdout, "%s\t%d\t%d\t%d\t%d\t%d\t%d", r.tn, r.tl, r.ts, r.te, r.ml, r.bl, r.mq);
      fprintf(stdout, "\t%s\n", r.cigar);
      */
      q = r.rev ? r.qe - 1 : r.qs;
      bin = kh_get(faHash, readmap, r.qn);
      if(bin == kh_end(readmap)) {
        fprintf(stderr, "ERROR: a read name in the PAF is apparently missing from the read FASTA.\n");
        return 1;
      }
      qi = kh_val(readmap, bin);
      qseq = kv_A(reads, qi).s;
      t = r.ts;
      bin = kh_get(faHash, refmap, r.tn);
      if(bin == kh_end(refmap)) {
        fprintf(stderr, "ERROR: a ref name in the PAF is apparently missing from the ref FASTA.\n");
        return 1;
      }
      ti = kh_val(refmap, bin);
      tseq = kv_A(refs, ti).s;
      val = 0;
      for(i = 0; i < strlen(r.cigar); i++) {
        if(r.cigar[i] < 58) { // a digit
          val = val * 10 + (r.cigar[i] - 48);
        } else {
          //fprintf(stderr, "%c: %u\n", r.cigar[i], val);
          if(r.cigar[i] == 'M' || r.cigar[i] == 'X' || r.cigar[i] == '=') { // consumes q and t
            for(j = 0; j < val; j++) {
              //fprintf(stderr, "qi: %u, t: %u, st: %u, q: %u\n", qi, t, st, q);
              if(t >= st && t < en)
                matrix[qi][t-st] = (r.rev ? compl[qseq[q]] : norm[qseq[q]]);
              q = q + (r.rev ? -1 : 1);
              t++;
            }
          }
          else if(r.cigar[i] == 'I' || r.cigar[i] == 'S') { // consumes q (INSERTION)
            q = q + (r.rev ? -1*val : val); // just skip it, we aren't recording insertions
          }
          else if(r.cigar[i] == 'D' || r.cigar[i] == 'N') { // consumes t (DELETION)
            for(j = 0; j < val; j++) { // fill in deletions (1s)
              if(t >= st && t < en)
                matrix[qi][t-st] = 45; // "-"
              t++;
            }
          }
          val = 0;
        }
      }
      ret = paf_read(p, &r);
    }
    paf_close(p);
    fprintf(stderr, "Closing file '%s'\n", paf_file);

    for(i = 0; i < kv_size(reads); i++) {
      fprintf(stderr, "%d |", i);
      for(j = 1500; j < 1650; j++) {
        fprintf(stderr, "%c", matrix[i][j-st]);
      }
      fprintf(stderr, "|\n");
    }
  }

  // find subset of reads that align across the entire region
  uint32_t n_full_reads = 0;
  kvec_t(uint32_t) full_reads;
  kv_init(full_reads);
  for(i = 0; i < kv_size(reads)-1; i++) {
    if(matrix[i][0] != ' ' && matrix[i][en-st-1] != ' ') { // alignment must cover the ENTIRE target region
      n_full_reads++;
      kv_push(uint32_t, full_reads, i);
    }
  }
  fprintf(stderr, "%u reads map fully from %d to %d (out of %u)\n", kv_size(full_reads), st, en, kv_size(reads));

  // ---------- pairwise distance ----------
  int k;
  uint16_t **dist = malloc((kv_size(full_reads)-1) * sizeof(uint16_t*));
  /*
   *   12345  dists:
   * A -0000  [[0, 0, 0, 0],
   * B --000   [0, 0, 0],
   * C ---00   [0, 0],
   * D ----0   [0]]
   * E -----
   */
  fprintf(stderr, "Computing distance matrix...\n");
  for(i = 0; i < kv_size(full_reads)-1; i++) {
    dist[i] = calloc((kv_size(full_reads) - (i+1)), sizeof(uint16_t));
    for(j = i+1; j < kv_size(full_reads); j++) {
      for(k = 0; k < en-st; k++) {
        if(matrix[kv_A(full_reads, i)][k] != matrix[kv_A(full_reads, j)][k]) {
          dist[i][j-i-1]++;
        }
      }
      //fprintf(stderr, "%d -- %d: %u\n", i, j, dist[i][j-i-1]);
    }
  }
  fprintf(stderr, "Agglomerative clustering...\n");
  hierarchy clusters = agglomerate_u16(dist, kv_size(full_reads));
  uint32_t cutoff = elbow_cutoff(clusters);
  fprintf(stderr, "cutoff: %u\n", cutoff);


  // ---------- output clusters ----------
  uint32_t *cluster_idx = calloc((kv_size(full_reads) + kv_size(clusters)), sizeof(uint32_t));
  // where the first |full_reads| are read IDs, and |full_reads| -> |full_reads| + |clusters| are higher level cluster IDs
  uint32_t cid = 1; // incremental cluster ID
  for(i = kv_size(clusters)-1; i >= 0; i--) {
  //for(i = 0; i < kv_size(clusters)-1; i++) {
    if(kv_size(clusters)-1-i <= cutoff) { // this merge should NOT be included, so each new cluster gets a new idx
    //if(i <= cutoff) { // this merge should NOT be included, so each new cluster gets a new idx
      fprintf(stderr, "keeping a(%u) <-> b(%u), dist %u: %f\n", kv_A(clusters, i).a, kv_A(clusters, i).b, i, kv_A(clusters, i).dist);
      cluster_idx[kv_A(clusters, i).a] = cluster_idx[i + kv_size(full_reads)];
      cluster_idx[kv_A(clusters, i).b] = cid++;
    } else { // these are subclusters, so they just inherit their idx
      cluster_idx[kv_A(clusters, i).a] = cluster_idx[i + kv_size(full_reads)];
      cluster_idx[kv_A(clusters, i).b] = cluster_idx[i + kv_size(full_reads)];
    }
  }

  fprintf(stderr, "%u total clusters\n", cid);

  for(i = 0; i < kv_size(full_reads); i++) {
    fprintf(stdout,  "%s\t%u\n", kv_A(readnames, kv_A(full_reads, i)).s, cluster_idx[i]);
  }

  // free: reads, readnames, readmap, refs, refmap, dists, matrix
  free(cluster_idx);
}
