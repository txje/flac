#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
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

#define NALLELES 7

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
  printf("Commands:\n");
  printf("  cluster: agglomerative clustering with elbow cutoff\n");
  printf("  consensus: generate a simple consensus from the given PAF\n");
  printf("  2snv: explicit denoising by finding linked variants with frequency > expected via binomial distribution\n");
  printf("Options:\n");
  printf("  -q: FASTA/Q[.gz] file with reads\n");
  printf("  -r: Reference FASTA/Q[.gz] or precomputed index file\n");
  printf("  -p, --paf: PAF alignment file, or '-' for stdin\n");
  printf("  -d, --distance-matrix: Write (or read, if it exists) distance matrix file\n");
  printf("  -b, --begin: begin position to evaluate (in single ref seq)\n");
  printf("  -e, --end: end position to evaluate (in single ref seq)\n");
  printf("\n");
  printf("  Alignment parameters:\n");
  printf("    -t: Threads (default: 1)\n");
  printf("    -s: Sequence type preset\n");
  printf("        map-ont: Oxford Nanopore (default)\n");
  printf("        map-pb:  Pacbio\n");
  printf("    -f, --align-fraction: Portion of a read that must align properly (defaults to --align-length threshold)\n");
  printf("    -l, --align-length: Minimum aligned bp (default: 100)\n");
  printf("    -a, --align-accuracy: Minimum accuracy of aligned portion of a read (default: 0.6)\n");
  printf("\n");
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
  { "distance-matrix",        required_argument, 0, 'd' },
  { "begin",                  required_argument, 0, 'b' },
  { "end",                    required_argument, 0, 'e' },
  { 0, 0, 0, 0}
};

// a shortcut to prevent numerical overflow when doing factorials
double choose(int n, int k, double a) {
  int i;
  for(i = 1; i <= k; i++) {
    a = (a * (n - (k - i))) / i;
  }
  return a;
}

// 2SNV uses exactly this kind of binomial CDF plus bonferroni correction
// but computes the probability of a specific variant explicitly as p = (#SNP1 * #SNP2) / (#NEITHER * #READS)
// so that the probability 1-binom_cdf(#SNP1&SNP2 - 1, #READS, p) must be <= 0.01 / choose(#LOCI, 2)

// chance of observing *up to* k successes among n trials with probability of success p
/*
double binom_cdf(int k, int n, float p) {
  double a = 0;
  int i;
  for(i = 0; i <= k; i++) {
    fprintf(stderr, "    i: %d, a: %f\n", i, a);
    fprintf(stderr, "      choose(%u, %u): %u\n", n, i, choose(n,i));
    fprintf(stderr, "      pow(%f, %u): %f\n", p, i, pow(p,i));
    fprintf(stderr, "      pow(1-%f, %u-%u): %f\n", p, n, i, pow(1-p,n-i));
    a = a + choose(n, i, pow(p, i) * pow(1-p, n-i));
  }
  return a;
}
*/

/*
 * Build a consensus sequence from the given matrix with l *columns* (nt)
 * - using only rows/cols where cluster_idx[i] == idx
 */
uint8_t* consensus(uint8_t **matrix, uint32_t l, uint32_t *full_reads, uint32_t n_full_reads, uint16_t *cluster_idx, uint32_t idx) {
  uint32_t i, j;
  uint32_t **alleles = malloc(l * sizeof(uint32_t*));
  for(i = 0; i < l; i++) {
    alleles[i] = calloc(7, sizeof(uint32_t)); // [unaligned], A, C, G, T, [del], N
  }
  for(i = 0; i < n_full_reads; i++) {
    if(cluster_idx == NULL || cluster_idx[i] == idx) {
      for(j = 0; j < l; j++) {
        alleles[j][matrix[full_reads[i]][j]]++;
      }
    }
  }
  uint8_t *cons = malloc(l * sizeof(uint8_t));
  uint8_t max;
  for(i = 0; i < l; i++) {
    max = 0;
    for(j = 1; j < 6; j++) {
      if(alleles[i][j] > alleles[i][max]) {
        max = j;
      }
    }
    cons[i] = max;
  }
  for(i = 0; i < l; i++) {
    free(alleles[i]);
  }
  free(alleles);
  return cons;
}

// convert a byte array (like the matrix or consensus) to a nucleotide string
char* c2str(uint8_t* c, uint32_t l) {
  char* s = malloc(l * sizeof(char));
  int i;
  for(i = 0; i < l; i++) {
    s[i] = " ACGT-N"[c[i]];
  }
  return s;
}


int main(int argc, char *argv[]) {
  mm_idxopt_t iopt;
  mm_mapopt_t mopt;

  char* read_fasta = NULL;
  char* ref_fasta = NULL;
  char* paf_file = NULL;
  char* preset = "map-ont";
  char* distance_matrix = NULL;
  float align_fraction = -1;
  float align_accuracy = 0.6;
  int align_length = 100;
  int n_threads = 1;
  int verbose = 0;
  int st = 0;
  int en = 0;

  /*
   * Probability of any given change where indices are:
   * unaligned = 0, A-T = 1-4, del = 5, N = 6
   * so, any simple SNP has 1% chance
   * any deletion has 4% chance
   * there should be NO unaligned (0) or Ns (6), so any that appear will be infinitely (?) unlikely to occur by chance
   */
  double PAIR_PROBS[] = 
    {0,    0,    0,    0,    0,    0,    0,
     0, 0.99, 0.01, 0.01, 0.01, 0.04,    0,
     0, 0.01, 0.99, 0.01, 0.01, 0.04,    0,
     0, 0.01, 0.01, 0.99, 0.01, 0.04,    0,
     0, 0.01, 0.01, 0.01, 0.99, 0.04,    0,
     0, 0.04, 0.04, 0.04, 0.04, 0.96,    0,
     0,    0,    0,    0,    0,    0,    0};

  // ---------- options ----------
  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "q:r:t:p:d:b:e:s:f:l:a:vh", long_options, &long_idx)) != -1) {
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
      case 'd':
        distance_matrix = optarg;
        break;
      case 'b':
        st = atoi(optarg);
        break;
      case 'e':
        en = atoi(optarg);
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      case '?':
        if (optopt == 'q' || optopt == 'r' || optopt == 't' || optopt == 'f' || optopt == 'a' || optopt == 'p' || optopt == 'l' || optopt == 's' || optopt == 'd' || optopt == 'b' || optopt == 'e')
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
        else if (long_idx == 7) distance_matrix = optarg; // --distance-matrix
        else if (long_idx == 8) st = atoi(optarg); // --begin
        else if (long_idx == 9) en = atoi(optarg); // --end
        break;
      default:
        usage();
        return 1;
    }
  }

  int index;
  char* command = NULL;
  for (index = optind; index < argc; index++) {
    if(index == optind) {
      command = argv[index];
    }
  }
  if(command == NULL) {
    usage();
    return 1;
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
  int i, j, k, l;

  // capitalize standard nucleotides, turn anything else into N
  /*
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
  */
  uint8_t norm[256];
  for(i = 0; i < 256; i++) {
    norm[i] = 6;
  }
  norm[65] = 1;
  norm[67] = 2;
  norm[71] = 3;
  norm[84] = 4;
  norm[97] = 1;
  norm[99] = 2;
  norm[103] = 3;
  norm[116] = 4;

  uint8_t compl[256];
  for(i = 0; i < 256; i++) {
    compl[i] = 6;
  }
  compl[65] = 4;
  compl[67] = 3;
  compl[71] = 2;
  compl[84] = 1;
  compl[97] = 4;
  compl[99] = 3;
  compl[103] = 2;
  compl[116] = 1;

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
  if(st == 0 && en == 0) { // GUESS?
    // for aav:
    st = 2250;
    en = 4450;
    // for simulated data:
    /*
    st=1100;
    en=2900;
    */
  }
  l = en-st;

  /*
   * matrix of reads x pileup, where uint8_t values are characters in this set:
   * unaligned = 0, A-T = 1-4, del = 5, N = 6
   */
  uint8_t **matrix = malloc(sizeof(uint8_t*) * kv_size(reads));
  for(i = 0; i < kv_size(reads); i++) {
    matrix[i] = calloc(l, sizeof(uint8_t)); // unaligned will be 0
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
                matrix[qi][t-st] = 5;
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

    if(verbose) {
      for(i = 0; i < kv_size(reads); i++) {
        fprintf(stderr, "%d |", i);
        for(j = st+l/2-75; j < st+l/2+75; j++) {
          fprintf(stderr, "%c", " ACGT-N"[matrix[i][j-st]]);
        }
        fprintf(stderr, "|\n");
      }
    }
  }

  // find subset of reads that align across the entire region
  uint32_t n_full_reads = 0;
  kvec_t(uint32_t) full_reads;
  kv_init(full_reads);
  for(i = 0; i < kv_size(reads)-1; i++) {
    if(matrix[i][0] && matrix[i][l-1]) { // alignment must cover the ENTIRE target region
      n_full_reads++;
      kv_push(uint32_t, full_reads, i);
    }
  }
  fprintf(stderr, "%u reads map fully from %d to %d (out of %u)\n", n_full_reads, st, en, kv_size(reads));
  
  if(strcmp(command, "cluster") == 0) {
    // ---------- pairwise distance ----------
    uint16_t **dist = malloc((n_full_reads-1) * sizeof(uint16_t*));
    /*
     *   12345  dists:
     * A -0000  [[0, 0, 0, 0],
     * B --000   [0, 0, 0],
     * C ---00   [0, 0],
     * D ----0   [0]]
     * E -----
     */

    FILE *df;
    uint32_t dm_size; // the squareform size (even though the distance matrix includes only n*(n-1)/2)
    uint8_t compute_dm = 0;
    size_t row_size;
    if(distance_matrix != NULL) {
      df = fopen(distance_matrix, "rb");
      if(!df || !fread(&dm_size, 4, 1, df)) {
        fprintf(stderr, "Didn't find anything in '%s'.\n", distance_matrix);
        compute_dm = 1;
      } else {
        if(dm_size != n_full_reads) {
          fprintf(stderr, "Number of full reads (%u) does not match the matrix size (%u)! We'll quit now to avoid overwriting the distance matrix.\n", n_full_reads, dm_size);
          return 1;
        }
        fprintf(stderr, "Reading distance matrix from '%s'...\n", distance_matrix);
        for(i = 0; i < dm_size-1; i++) {
          row_size = (dm_size - (i+1)) * sizeof(uint16_t);
          dist[i] = malloc(row_size);
          fread(dist[i], row_size, 1, df);
        }
      }
      if(df)
        fclose(df);
    } else {
      compute_dm = 1;
    }

    if(compute_dm) {
      dm_size = n_full_reads;
      fprintf(stderr, "Computing distance matrix...\n");
      for(i = 0; i < dm_size-1; i++) {
        dist[i] = calloc((dm_size - (i+1)), sizeof(uint16_t));
        for(j = i+1; j < dm_size; j++) {
          for(k = 0; k < l; k++) {
            if(matrix[kv_A(full_reads, i)][k] != matrix[kv_A(full_reads, j)][k]) {
              dist[i][j-i-1]++;
            }
          }
          //fprintf(stderr, "%d -- %d: %u\n", i, j, dist[i][j-i-1]);
        }
      }
      if(distance_matrix != NULL) {
        fprintf(stderr, "Writing distance matrix (%u reads) to '%s'\n", dm_size, distance_matrix);
        df = fopen(distance_matrix, "wb");
        if(!df) {
          fprintf(stderr, "Failed to open '%s' for writing!\n", distance_matrix);
          return 1;
        }
        fwrite(&dm_size, 4, 1, df); // write 4-byte matrix size
        for(i = 0; i < dm_size-1; i++) {
          row_size = (dm_size - (i+1)) * sizeof(uint16_t);
          fwrite(dist[i], row_size, 1, df);
        }
        fclose(df);
      }
    }
    fprintf(stderr, "Agglomerative clustering...\n");
    hierarchy clusters = agglomerate_u16(dist, n_full_reads);
    uint32_t cutoff = elbow_cutoff(clusters);
    fprintf(stderr, "cutoff: %u\n", cutoff);


    // ---------- output clusters ----------
    uint32_t *cluster_idx = calloc((n_full_reads + kv_size(clusters)), sizeof(uint32_t));
    // where the first |full_reads| are read IDs, and |full_reads| -> |full_reads| + |clusters| are higher level cluster IDs
    uint32_t cid = 1; // incremental cluster ID
    for(i = kv_size(clusters)-1; i >= 0; i--) {
      if(kv_size(clusters)-1-i <= cutoff) { // this merge should NOT be included, so each new cluster gets a new idx
        if(verbose) {
          fprintf(stderr, "Splitting bicluster %u: a(%u) <-> b(%u), dist: %f\n", i, kv_A(clusters, i).a, kv_A(clusters, i).b, kv_A(clusters, i).dist);
        }
        cluster_idx[kv_A(clusters, i).a] = cluster_idx[i + n_full_reads];
        cluster_idx[kv_A(clusters, i).b] = cid++;
      } else { // these are subclusters, so they just inherit their idx
        cluster_idx[kv_A(clusters, i).a] = cluster_idx[i + n_full_reads];
        cluster_idx[kv_A(clusters, i).b] = cluster_idx[i + n_full_reads];
      }
    }

    fprintf(stderr, "%u total clusters\n", cid);

    for(i = 0; i < n_full_reads; i++) {
      fprintf(stdout,  "%s\t%u\n", kv_A(readnames, kv_A(full_reads, i)).s, cluster_idx[i]);
    }

    // ---------- clean up clustering memory ----------
    for(i = 0; i < n_full_reads-1; i++) {
      free(dist[i]);
    }
    free(dist);
    free(cluster_idx);

  } else if (strcmp(command, "2snv") == 0) {
    // ---------- variant denoising a la 2SNV ----------
    uint16_t *cluster_idx = calloc(n_full_reads, sizeof(uint16_t)); // all reads start in cluster 0

    kvec_t(uint8_t*) cluster_seqs; // vector of cluster consensus sequences
    kv_init(cluster_seqs);
    kv_push(uint8_t*, cluster_seqs, consensus(matrix, l, full_reads.a, n_full_reads, cluster_idx, 0));

    uint16_t cid = 1; // incremental cluster index
    uint8_t keep_splitting = 1; // a flag to mark when we're done iterating
    uint32_t one, two; // indices of the two (possibly linked) loci
    uint32_t best_1, best_2, best_c, best_ct;
    double prob, p;
    uint32_t ***alleles; // for each cluster, a 6x6 matrix of allele-pair counts
    uint8_t a, b, major_a, major_b, minor_a, minor_b; // alleles
    //uint32_t *a_ct = malloc(NALLELES * sizeof(uint32_t)); // allele subset counts
    //uint32_t *b_ct = malloc(NALLELES * sizeof(uint32_t));
    uint32_t *cluster_sizes;
    uint32_t n_loci_pairs = (uint32_t)choose(l, 2, 1);
    uint32_t majority_ct; // the # of reads with the majority alleles in a single cluster, single pair of loci
    // successively split
    while(keep_splitting) {
      // reset allele counts
      alleles = malloc(cid * sizeof(uint32_t**));
      for(i = 0; i < cid; i++) {
        alleles[i] = malloc(NALLELES * sizeof(uint32_t*));
        for(j = 0; j < NALLELES; j++) {
          alleles[i][j] = malloc(NALLELES * sizeof(uint32_t));
        }
      }
      fprintf(stderr, "%u clusters\n", cid);

      // recompute cluster sizes
      cluster_sizes = calloc(cid, sizeof(uint32_t));
      for(i = 0; i < n_full_reads; i++) {
        cluster_sizes[cluster_idx[i]]++;
      }

      // precompute significant read cutoffs per cluster per allele-pair (allele1 a->b (7x7) x allele2 a->b (7x7) * #clusters), so 49^2 * clusters
      uint32_t *cutoffs = malloc((cid*NALLELES*NALLELES*NALLELES*NALLELES) * sizeof(uint32_t));
      for(i = 0; i < cid; i++) {
        fprintf(stderr, "cluster %u: %u reads\n", i, cluster_sizes[i]);
        for(a = 0; a < NALLELES*NALLELES; a++) { // represents a single variant in the space of 7x7 changes
          for(b = 0; b < NALLELES*NALLELES; b++) { // represents a single variant in the space of 7x7 changes
            p = 0;
            for(j = 0; j < cluster_sizes[i]; j++) {
              //fprintf(stderr, "exactly %u: %.20f\n", j, (choose(cluster_sizes[i], j, pow(pair_prob, j) * pow(1-pair_prob, cluster_sizes[i]-j))));
              p = p + (choose(cluster_sizes[i], j, pow(PAIR_PROBS[a]*PAIR_PROBS[b], j) * pow(1-(PAIR_PROBS[a]*PAIR_PROBS[b]), cluster_sizes[i]-j)));
              //fprintf(stderr, "binom_cdf(%u, %u, %f) = %.10f\n", j, cluster_sizes[i], PAIR_PROBS, p);
              if(1-p < 0.01 / n_loci_pairs) {
                cutoffs[i*NALLELES*NALLELES*NALLELES*NALLELES+a*NALLELES*NALLELES+b] = j;
                break;
              }
            }
            fprintf(stderr, "  cutoff for %u -> %u AND %u -> %u: %u\n", a/NALLELES, a%NALLELES, b/NALLELES, b%NALLELES, cutoffs[i*NALLELES*NALLELES*NALLELES*NALLELES+a*NALLELES*NALLELES+b]);
          }
        }
      }

      // check if any cluster (idx < cid) has 2 significantly linked SNVs
      best_ct = 0; // best total number of linked reads
      best_1 = 0;
      best_2 = 1;
      best_c = 0;
      for(one = 0; one < l-1; one++) {
        //fprintf(stderr, "checking %u vs. all\n", one);
        for(two = one+1; two < l; two++) {
          // reset allele counts
          for(i = 0; i < cid; i++) {
            for(j = 0; j < NALLELES; j++) {
              for(k = 0; k < NALLELES; k++) {
                alleles[i][j][k] = 0;
              }
            }
          }

          //fprintf(stderr, "checking loci %u and %u\n", one, two);
          // build allele matrix for all cids at once
          for(i = 0; i < n_full_reads; i++) {
            a = matrix[kv_A(full_reads,i)][one];
            b = matrix[kv_A(full_reads,i)][two];
            alleles[cluster_idx[i]][a][b]++;
          }
          // test each pair for each cid
          for(i = 0; i < cid; i++) {
            //fprintf(stderr, "  checking cluster %u\n", i);
            // compute per-row and per-col totals
            /*
            for(a = 0; a < NALLELES; a++) {
              a_ct[a] = 0;
              b_ct[a] = 0;
            }
            for(a = 0; a < NALLELES; a++) {
              for(b = 0; b < NALLELES; b++) {
                a_ct[a] += alleles[i][a][b];
                b_ct[b] += alleles[i][a][b];
              }
            }
            */
            /*
            fprintf(stderr, "cluster has %u reads\n", cluster_size);
            for(a = 0; a < NALLELES; a++) {
              for(b = 0; b < NALLELES; b++) {
                fprintf(stderr, "  %u/%u: %u\n", a, b, alleles[i][a][b]);
              }
            }
            */
            majority_ct = 0;
            for(a = 0; a < NALLELES; a++) {
              for(b = 0; b < NALLELES; b++) {
                if(alleles[i][a][b] == 0) continue;
                if(alleles[i][a][b] > majority_ct) {
                  majority_ct = alleles[i][a][b];
                  major_a = a;
                  major_b = b;
                }
              }
            }
            for(a = 0; a < NALLELES; a++) {
              for(b = 0; b < NALLELES; b++) {
                if(a == major_a || b == major_b || alleles[i][a][b] == 0) continue;
                // 2SNV computes the probability of a specific variant explicitly as p = (#SNP1 * #SNP2) / (#NEITHER * #READS)
                // so that the probability of observing is 1-binom_cdf(#SNP1&SNP2 - 1, #READS, p) and must be <= 0.01 / choose(#LOCI, 2)
                // this shit doesn't work --
                /*
                if(a_ct[a] > alleles[i][a][b]) {
                  if(b_ct[b] > alleles[i][a][b]) {
                    prob = ((a_ct[a]-alleles[i][a][b]) * (b_ct[b]-alleles[i][a][b])) / (double)((cluster_size-a_ct[a]-b_ct[b]+alleles[i][a][b]) * cluster_size);
                  } else {
                    prob = ((a_ct[a]-alleles[i][a][b])) / (double)((cluster_size-a_ct[a]-b_ct[b]+alleles[i][a][b]) * cluster_size);
                  }
                } else {
                  prob = ((b_ct[b]-alleles[i][a][b])) / (double)((cluster_size-a_ct[a]-b_ct[b]+alleles[i][a][b]) * cluster_size);
                }
                */
                /*
                prob = 0.01;
                //fprintf(stderr, "stat check for %u reads of %u with prob = %f\n", alleles[i][a][b], cluster_size, prob);
                //TODO: chance uses ALL full reads, not the current cluster_size, so it's WRONG
                p = 1 - chance[alleles[i][a][b]-1]; //binom_cdf(alleles[i][a][b] - 1, cluster_size, prob, chose);
                //fprintf(stderr, "p: %f\n", p);
                if(p < best_p) {
                  fprintf(stderr, "New best p-val: %f at cluster %u, %u : %u\n", p, i, one, two);
                  best_p = p;
                  best_1 = one;
                  best_2 = two;
                  best_c = cid;
                }
                */
                if(alleles[i][a][b] >= cutoffs[i*NALLELES*NALLELES*NALLELES*NALLELES+a*NALLELES*NALLELES*NALLELES+major_a*NALLELES*NALLELES+b*NALLELES+major_b] && alleles[i][a][b] > best_ct) {
                  fprintf(stderr, "  %u -> %u AND %u -> %u has %u reads, needs %u\n", a, major_a, b, major_b, alleles[i][a][b], cutoffs[i*NALLELES*NALLELES*NALLELES*NALLELES+a*NALLELES*NALLELES*NALLELES+major_a*NALLELES*NALLELES+b*NALLELES+major_b]);
                  best_ct = alleles[i][a][b];
                  best_1 = one;
                  best_2 = two;
                  best_c = i;
                  minor_a = a;
                  minor_b = b;
                }
              }
            }
          }
        }
      }
      if(best_ct == 0) {
        fprintf(stderr, "No more linked SNVs, finishing up...\n");
        keep_splitting = 0;
      } else {
        fprintf(stderr, "Best pair: %u (%u) and %u (%u) in cluster %u: reads = %u\n", best_1, minor_a, best_2, minor_b, best_c, best_ct);
        // split best_c by this pair
        // put every read that has the BOTH new alleles into a new cluster
        for(i = 0; i < n_full_reads; i++) {
          if(cluster_idx[i] == best_c) {
            if(matrix[kv_A(full_reads,i)][best_1] == minor_a && matrix[kv_A(full_reads,i)][best_2] == minor_b) {
              cluster_idx[i] = cid;
            }
          }
        }
        // make consensus for new cluster
        kv_push(uint8_t*, cluster_seqs, consensus(matrix, l, full_reads.a, n_full_reads, cluster_idx, cid));
        // and remake consensus for split cluster
        free(kv_A(cluster_seqs, best_c));
        kv_A(cluster_seqs, best_c) = consensus(matrix, l, full_reads.a, n_full_reads, cluster_idx, best_c);
        // reassign ALL reads to its most similar cluster
        uint32_t diff, best_diff;
        for(i = 0; i < n_full_reads; i++) {
          best_c = 0; // REPURPOSING best_c here
          best_diff = l;
          for(j = 0; j <= cid; j++) {
            diff = 0;
            for(k = 0; k < l; k++) {
              if(matrix[kv_A(full_reads, i)][k] != kv_A(cluster_seqs, j)[k]) {
                diff++;
              }
            }
            if(diff < best_diff) {
              best_c = j;
              best_diff = diff;
            }
          }
          cluster_idx[i] = best_c;
        }
        // recompute consensus sequences for all clusters
        for(j = 0; j <= cid; j++) {
          free(kv_A(cluster_seqs, j));
          kv_A(cluster_seqs, j) = consensus(matrix, l, full_reads.a, n_full_reads, cluster_idx, j);
        }
      }

      // clean up temporary splitting data for next round
      for(i = 0; i < cid; i++) {
        for(j = 0; j < NALLELES; j++)
          free(alleles[i][j]);
        free(alleles[i]);
      }
      free(alleles);
      free(cluster_sizes);
      free(cutoffs);

      if(keep_splitting)
        cid++;
    }
    // print then free cluster_seqs
    //
    // cluster info is in cluster_idx[full_read_idx] and cluster_seqs(vec)
    fprintf(stderr, "%u clusters found\n", cid);
    for(i = 0; i < n_full_reads; i++) {
      fprintf(stdout,  "%s\t%u\n", kv_A(readnames, kv_A(full_reads, i)).s, cluster_idx[i]);
    }
    free(cluster_idx);

  } else if (strcmp(command, "consensus") == 0) {
    // ---------- generate a simple consensus using all full-length reads ----------
    uint8_t *cons = consensus(matrix, l, full_reads.a, n_full_reads, NULL, 0);
    char *conseq = c2str(cons, l);
    free(cons);
    fprintf(stdout, ">consensus\n");
    for(j = 0; j < en-st; j++) {
      fprintf(stdout, "%c", conseq[j]);
      if((j+1) % 80 == 0) {
        fprintf(stdout, "\n");
      }
    }
    free(conseq);
  } else {
    fprintf(stderr, "Command '%s' not recognized.\n", command);
    return 1;
  }


  // ---------- clean up general memory ----------
  for(i = 0; i < kv_size(reads); i++) {
    free(matrix[i]);
  }
  free(matrix);

  for(bin = kh_begin(refmap); bin != kh_end(refmap); bin++) {
    if(kh_exist(refmap, bin))
      free((char*)kh_key(refmap, bin));
  }
  kh_destroy(faHash, refmap);
  for(i = 0; i < kv_size(refs); i++) {
    free(kv_A(refs, i).s);
  }
  kv_destroy(refs);
  for(i = 0; i < kv_size(reads); i++) {
    free(kv_A(reads, i).s);
    free(kv_A(readnames, i).s);
  }
  kv_destroy(reads);
  kv_destroy(readnames);
  kh_destroy(faHash, readmap); // name keys were freed from readnames vector
  kv_destroy(full_reads); // contains only ints

  return 0;
}
