/*
 * https://github.com/lh3/miniasm/blob/master/paf.h
 * commit ff3f28a
 */

#ifndef PAF_PAF_H
#define PAF_PAF_H

#include <stdint.h>
#include <sys/types.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
  size_t l, m;
  char *s;
} kstring_t;
#endif

typedef struct {
  void *fp;
  kstring_t buf;
} paf_file_t;

typedef struct {
  const char *qn, *tn; // these point to the input string; NOT allocated
  uint32_t ql, qs, qe, tl, ts, te;
  uint32_t ml:31, rev:1, bl;
  uint8_t mq;
  uint8_t *cigar; // length == bl, uncompressed CIGAR string giving edits where {0:match, 1:ins, 2:del}, no soft/hard-clipping
} paf_rec_t;

#ifdef __cplusplus
extern "C" {
#endif

paf_file_t *paf_open(const char *fn);
int paf_close(paf_file_t *pf);
int paf_read(paf_file_t *pf, paf_rec_t *r);

#ifdef __cplusplus
}
#endif

#endif
