#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "incl/klib/kvec.h"
#include "hierarch.h"

/*
 * Lance-Williams distance calculation a la Wikipedia
 * *assumes* squared euclidean distance, which we do not have...
 * clusters a and b will be merged, so compute distance (a U b) to all other clusters (c)
 */
int compute_dist(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b) {
  uint32_t c;
  float dist;
  for(c = 0; c < n_clust; c++) {
    if(c != a && c != b && sizes[c] > 0) {
      float tot =  (float)(sizes[a] + sizes[b] + sizes[c]);
      dist = ((sizes[a] + sizes[c]) / tot * dists[a][c]) \
           + ((sizes[b] + sizes[c]) / tot * dists[b][c]) \
           - ((sizes[c]) / tot * dists[a][b]);
      //fprintf(stderr, "%u/%u : %u dist = %f\n", a, b, c, dist);
      dists[n_clust][c] = dist;
      dists[c][n_clust] = dist;
    }
  }
}

/*
 * average pairwise distance
 */
int compute_dist_pairwise(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b) {
  uint32_t c;
  float dist;
  for(c = 0; c < n_clust; c++) {
    if(c != a && c != b && sizes[c] > 0) {
      float dist;
      // do this somehow differently...
      fprintf(stderr, "%u/%u : %u dist = %f\n", a, b, c, dist);
      dists[n_clust][c] = dist;
      dists[c][n_clust] = dist;
    }
  }
}

int find_closest(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t *a, uint32_t *b) {
  uint32_t i, j;
  uint32_t min = 123456789;
  for(i = 0; i < n_clust; i++) {
    if(sizes[i] == 0) continue;
    for(j = 0; j < n_clust; j++) {
      if(i == j || sizes[j] == 0) continue;
      if(min == 123456789 || dists[i][j] < min) {
        *a = i;
        *b = j;
        min = dists[i][j];
      }
    }
  }
}

int merge(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b) {
  compute_dist(dists, sizes, n_clust, a, b);
  sizes[n_clust] = sizes[a] + sizes[b];
  sizes[a] = 0;
  sizes[b] = 0;
}

void print_square_matrix(FILE *o, float **a, uint32_t n) {
  int i, j;
  for(j = 0; j < n; j++) {
    fprintf(o, "--");
  }
  fprintf(o, "\n");
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      fprintf(o, "%.0f ", a[i][j]);
    }
    fprintf(o, "\n");
  }
  for(j = 0; j < n; j++) {
    fprintf(o, "--");
  }
  fprintf(o, "\n");
}

hierarchy ag(float **cluster_dists, uint32_t *sizes, uint32_t n_items) {
  uint32_t i;

  // iteratively merge clusters
  hierarchy clusters;
  kv_init(clusters);
  uint32_t n_clust = n_items;
  uint32_t a, b;
  while(n_clust < n_items * 2 - 1) {
    //fprintf(stderr, "dists:\n");
    //print_square_matrix(stderr, cluster_dists, n_clust);
    find_closest(cluster_dists, sizes, n_clust, &a, &b);

    bicluster clust;
    clust.a = a;
    clust.b = b;
    clust.dist = cluster_dists[a][b];
    kv_push(bicluster, clusters, clust);

    //fprintf(stdout, "merging clusters %u and %u with distance %f\n", a, b, cluster_dists[a][b]);
    merge(cluster_dists, sizes, n_clust, a, b);
    n_clust++;
  }

  return clusters;
}

// the value this returns means the break is between biclusters idx (n - 1 - <elbow_cutoff>) and (n - 2 - <elbow_cutoff>)
int elbow_cutoff(hierarchy clusters) {
  uint32_t i;

  // compute elbow cutoff
  float *deriv = malloc(kv_size(clusters) * sizeof(float));
  for(i = 0; i < kv_size(clusters)-1; i++) {
    deriv[i] = kv_A(clusters, kv_size(clusters)-1-i).dist - kv_A(clusters, kv_size(clusters)-2-i).dist; // discrete derivative
    //fprintf(stderr, "%.2f ", deriv[i]);
    if(i > 0) {
      deriv[i-1] = deriv[i] - deriv[i-1]; // and second derivative trailing behind
      //fprintf(stderr, "%.2f ", deriv[i-1]);
    }
  }
  uint32_t argmin = 0;
  for(i = 1; i < kv_size(clusters)-2; i++) {
    if(deriv[i] < deriv[argmin]) argmin = i;
  }

  free(deriv);
  return argmin;
}

hierarchy agglomerate_u16(uint16_t **dists, uint32_t n_items) {
  // initialize cluster sizes and distances
  uint32_t *sizes = malloc(n_items * 2 * sizeof(uint32_t));
  float **cluster_dists = malloc(n_items * 2 * sizeof(float*));
  uint32_t i, j;
  for(i = 0; i < n_items * 2; i++) {
    sizes[i] = 1;
    cluster_dists[i] = calloc(n_items * 2, sizeof(float));
  }
  /* if a full symmetrical matrix:
  for(i = 0; i < n_items; i++) {
    for(j = 0; j < n_items; j++) {
      cluster_dists[i][j] = dists[i][j];
    }
  }
  */

  // if a partial matrix (each row is successively smaller)
  /*
   *   12345  dists:
   * A -0000  [[0, 0, 0, 0],
   * B --000   [0, 0, 0],
   * C ---00   [0, 0],
   * D ----0   [0]]
   * E -----
   */
  fprintf(stderr, "Expanding distance matrix...\n");
  for(i = 0; i < n_items; i++) {
    for(j = 0; j < n_items - i - 1; j++) {
      cluster_dists[i][j+i+1] = (float)dists[i][j];
      cluster_dists[j+i+1][i] = (float)dists[i][j];
    }
  }
  fprintf(stderr, "Aggregating...\n");
  hierarchy res = ag(cluster_dists, sizes, n_items);

  free(cluster_dists);
  free(sizes);
  return res;
}

hierarchy agglomerate(float **dists, uint32_t n_items) {

  // initialize cluster sizes and distances
  uint32_t *sizes = malloc(n_items * 2 * sizeof(uint32_t));
  float **cluster_dists = malloc(n_items * 2 * sizeof(uint32_t*));
  uint32_t i, j;
  for(i = 0; i < n_items * 2; i++) {
    sizes[i] = 1;
    cluster_dists[i] = calloc(n_items * 2, sizeof(uint32_t));
  }
  /* if a full symmetrical matrix:
  for(i = 0; i < n_items; i++) {
    for(j = 0; j < n_items; j++) {
      cluster_dists[i][j] = dists[i][j];
    }
  }
  */
  // if a partial matrix (each row is successively smaller)
  for(i = 0; i < n_items; i++) {
    for(j = 0; j < n_items - i - 1; j++) {
      cluster_dists[i][j+i+1] = dists[i][j];
      cluster_dists[j+i+1][i] = dists[i][j];
    }
  }

  hierarchy res = ag(cluster_dists, sizes, n_items);

  free(cluster_dists);
  free(sizes);
  return res;
}

