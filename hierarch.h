typedef struct bicluster {
  uint32_t a;
  uint32_t b;
  float dist;
} bicluster;

typedef kvec_t(bicluster) hierarchy;

int compute_dist(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b);
int find_closest(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t *a, uint32_t *b);
int merge(float **dists, uint32_t *sizes, uint32_t n_clust, uint32_t a, uint32_t b);
hierarchy agglomerate(float **dists, uint32_t n_items);
hierarchy agglomerate_u16(uint16_t **dists, uint32_t n_items);
int elbow_cutoff(hierarchy clusters);
