typedef struct bicluster {
  uint32_t a;
  uint32_t b;
  float dist;
} bicluster;

int compute_dist(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b);
int find_closest(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t *a, uint32_t *b);
int merge(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b);
int agglomerate(float **dists, uint32_t n_items);
