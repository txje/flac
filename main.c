#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "hierarch.h"

int main(int argc, char *argv[]) {
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
