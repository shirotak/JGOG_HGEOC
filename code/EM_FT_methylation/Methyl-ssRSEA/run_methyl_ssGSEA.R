# Custom script From Methyl-ssRSEA (PMID:34241544)
# $1 input tsv
# $2 region set bed (.gz is OK) 
# $3 probe coord file, 
# $4 output file path
# $5 score name, must
## check index name of probe ($1), default:'ID_REF', if needed, change Line 100

# input ------------------------------------------------------------------
args=commandArgs(trailingOnly=TRUE)
beta_value_path = args[1]
region_set_path = args[2]
probe_coord_path = args[3]
output_path = args[4]
score_name = args[5]

# packages --------------------------------------------------------------- 
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(fastmatch)
library(data.table)
# functions ---------------------------------------------------------------
mk_CpG_flag <- function(probe_coord, region_set){
  
  region4bind <- region_set %>%
    select(chr,start,end) %>%
    arrange(start) %>%
    mutate(id = paste0("rs_", row_number())) %>%
    tidyr::gather(key = start_end, value = pos, start, end) %>%
    mutate(start_end =
             if_else(start_end == "start", 1L, if_else(start_end == "end", -1L, as.integer(NA)))) %>%
    select(id, start_end, chr, pos) %>%
    arrange(chr, pos)
  
  CpG4bind <- probe_coord %>%
    mutate(start_end = 0) %>%
    select(id, start_end, chr, pos) %>%
    arrange(chr, pos)
  
  CpG_flag <- rbind(CpG4bind, region4bind) %>%
    arrange(chr, pos, start_end) %>%
    mutate(csum = cumsum(start_end)) %>%
    dplyr::filter(start_end == 0) %>%
    mutate(flag = pmin(csum, 1)) %>%
    select(id, flag) %>%
    arrange(id)
  
  return(CpG_flag)
}


mk_CpG_flag_dev01 <- function(probe_coord, region_set){
  bed <- region_set %>% mutate(flag = 1) %>% as.data.table
  CpG_flag <- probe_coord %>%
    as.data.table %>%
    bed[., on = .(chr, start < pos, end > pos), nomatch = 0, .(id, flag)] %>% # "nomatch=0" returns a inner join
    distinct(id, flag)
  return(CpG_flag)
}

calc_score <- function(R_vec, geneSets, allProbes, alpha){
  n <- length(R_vec)
  
  names(R_vec) = allProbes
  R_vec <- sort(R_vec, decreasing = TRUE)
  Ra_vec = abs(R_vec)^alpha
  # names(Ra_vec) = allProbes # 不要では？
  
  idxs <- sort(fmatch(unlist(geneSets), unlist(names(R_vec))))
  k <- length(idxs)
  stepCDFinGeneSet  <- sum(Ra_vec[idxs]*(n-idxs+1))/sum(Ra_vec[idxs])
  stepCDFoutGeneSet <- (n*(n+1)/2 - sum(n-idxs+1))/(n-k)
  
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
  
  return(walkStat)
}


# main --------------------------------------------------------------------
alpha = 1

print(paste("1. start", Sys.time()))
### input ###
beta_value <- read_tsv(beta_value_path)
region_set <- read_tsv(region_set_path, col_names=FALSE)[,1:3]
colnames(region_set) = c("chr", "start", "end")
probe_coord <- read_csv(probe_coord_path) %>% select(id, chr, pos)

print(paste("2. region set flag", Sys.time()))
### attach flag to CpG-probes in the region set ###
# CpG_flag <- mk_CpG_flag(probe_coord, region_set)
# geneSets <- CpG_flag %>% arrange(id) %>% dplyr::filter(flag == 1) %>% .$id %>% as.character() # list
geneSets <- mk_CpG_flag_dev01(probe_coord, region_set) %>% arrange(id) %>% dplyr::filter(flag == 1) %>% .$id %>% as.character() # list

print(paste("3: calc", Sys.time()))
### calculate scores ###
X <- beta_value %>% arrange(ID_REF) %>% column_to_rownames("ID_REF") %>% as.matrix()
# R <- apply(X, 2, function(x) as.integer(rank(x)))
# allProbes <- rownames(X)
# ssRSEA_scores <- apply(R, 2, calc_score, geneSets, allProbes, alpha)
ssRSEA_scores <- apply(X, 2, function(x){
  names(x) <- rownames(X)
  print(paste(sum(is.na(x)), sum(!is.na(x))))
  #print(head(x[is.na(x)]))
  x <- x[!is.na(x)] # remove NA
  x_allProbes <- names(x)
  x_geneSets <- geneSets[geneSets %in% x_allProbes]
  print(paste(length(x_geneSets), length(geneSets)))
  R_vec <- as.integer(rank(x))
  #print(head(x_allProbes))
  #print(head(x_geneSets))
  #print(head(R_vec))
  o <- calc_score(R_vec, x_geneSets, x_allProbes, alpha)
  return(o)
})

print(paste("4: output", Sys.time()))
### output ###
dfw=as.data.frame(ssRSEA_scores)
colnames(dfw)=score_name
write.table(dfw, output_path, quote = F,sep='\t',row.names = T,col.names = NA)

print(paste("5: end", Sys.time()))
