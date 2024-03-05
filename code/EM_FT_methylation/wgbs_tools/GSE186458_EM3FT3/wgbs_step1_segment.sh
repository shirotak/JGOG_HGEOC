#!/usr/bin/bash
# go to *beta directory

 ~/local/softwares/wgbs_tools/wgbstools segment --betas ./*beta --min_cpg 3 --max_bp 2000 -o blocks.small.bed

 ~/local/softwares/wgbs_tools/wgbstools index blocks.small.bed

