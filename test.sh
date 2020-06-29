#!/bin/bash

reads=test/aav.1000-3000.d0.01.c16x100.nanosim.fasta
ref=test/aav_1.ncbi.fasta
paf=test/ns_0.01x16_to_aav.paf

#./flac cluster -q $reads -r $ref -p $paf -d aav_nanosim.dist -b 1100 -e 2900 > aav_nanosim_cluster_test.tsv

echo
echo "Expect 16 clusters x 100 reads each"
echo "-----------------------------------"
awk '{
  if(!($2 in clusts)) {
    nclust++;
  }
  clusts[$2]++
} END {
  print nclust" clusters";
  for(c in clusts) {
    print "cluster "c": "clusts[c]" reads";
  }
}' aav_nanosim_cluster_test.tsv

./flac consensus -q $reads -r $ref -p $paf -b 1100 -e 2900 > aav_nanosim_consensus.fasta

#./flac 2snv -q $reads -r $ref -p $paf -b 1100 -e 2900 > aav_nanosim_2snv_test.tsv
