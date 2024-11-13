#!/usr/bin/bash
# $1 input block.small.bed.gz
# $2 group.csv
# $3 output dir, existing
# go to *beta directory

~/local/softwares/wgbs_tools/wgbstools find_markers --blocks_path $1 --groups_file $2 -o $3 --betas *beta --only_hyper
