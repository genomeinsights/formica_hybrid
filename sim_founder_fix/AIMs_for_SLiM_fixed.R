# Corrected version of SLiM/R_scripts/AIMs_for_SLiM.R
# (github.com/bcportinha/Replicate-hybrid-evolution)
#
# Output format is unchanged (one file per chromosome, tab-separated, no header):
#   chr_id  pos_0based  cluster  pm10_aq  pm10_pol
# so SpecIAnt_rufa_genome_LDclusters_demo_server.slim reads it as before.
#
# What changed (see NOTE_for_Beatriz.md):
#  1. BUG FIX: cluster means were computed as mean(aims_raw$faqu_fix[cluster_i]),
#     i.e. a WITHIN-chromosome row index applied to the GENOME-WIDE table. Every
#     chromosome except chromosome 1 therefore received chromosome 1's fixation
#     levels. Values are now looked up by marker ("chr:pos"), never by row index.
#  2. The trimming of AIMs past the chromosome end was computed (aims_trimmed)
#     but never used -- the untrimmed table was written. Trimming is now applied
#     to what is written.
#  3. Clusters are numbered per chromosome exactly as before (LD clusters first,
#     in order of first appearance in group_info, then independent SNPs), so the
#     SLiM side needs no change.

options(scipen=999)

setwd("./rufa_assembly/")

aims_raw <- read.table("species_diagnostic_marker_fixation_levels.tsv", header=T, sep="\t")

aims_raw$chr_id <- as.integer(gsub("chromosome_", "", aims_raw$chromosome))
aims_raw$marker <- paste(aims_raw$chr_id, ":", aims_raw$position, sep="")
stopifnot(!anyDuplicated(aims_raw$marker))

# 1-based -> 0-based, consistent with the recombination map conversion
aims_raw$pos_0based <- aims_raw$position - 1

# chromosome lengths - subtract 1 to make them 0-based
Ls = c(16449812, 18558483, 15805178, 16666556, 14957794, 10925039, 13372598, 13028660, 10584038, 13832892, 11830133,
       11260715, 11449430, 7730741, 10965718, 11699404, 11021207, 8644910, 7297899, 9820556, 10314430, 7671390,
       NA, 6269374, 13114028, 9969264, 7047176) - 1

# LD clusters from Petri
ld_clusters <- readRDS("group_info_new.rds")
ld_clusters$chr_id <- as.integer(gsub("chromosome_", "", ld_clusters$chromosome))
ld_clusters$marker <- paste(ld_clusters$chr_id, ":", ld_clusters$position, sep="")
stopifnot(!anyDuplicated(ld_clusters$marker))

# remove sites not present in the AIMs data
ld_clusters <- ld_clusters[ld_clusters$marker %in% aims_raw$marker, ]


for(chr in sort(unique(aims_raw$chr_id))){

  # this chromosome only, in physical order
  tmp <- aims_raw[aims_raw$chr_id == chr, ]
  tmp <- tmp[order(tmp$position), ]

  # are there any AIMs past the end of the chromosome?
  remove_i <- which(tmp$pos_0based > Ls[chr])
  if(length(remove_i) > 0){
    print(paste(chr, " has AIMs past its end! Trimming ", length(remove_i), " AIMs!", sep=""))
    tmp <- tmp[-remove_i, ]
  }

  ### identify sites in the same LD cluster
  ld_chr <- ld_clusters[ld_clusters$chr_id == chr, ]

  tmp$cluster <- NA
  tmp$aq_fix <- NA
  tmp$pol_fix <- NA

  cluster_ids <- unique(ld_chr$group_id)
  for(i in seq_along(cluster_ids)){

    # markers of this LD cluster -> rows of tmp (matched by marker, not by index)
    cluster_markers <- ld_chr$marker[ld_chr$group_id == cluster_ids[i]]
    cluster_i <- which(tmp$marker %in% cluster_markers)
    if(length(cluster_i) == 0) next

    tmp$cluster[cluster_i] <- i

    # mean fixation levels over THIS cluster's own markers
    tmp$aq_fix[cluster_i] <- mean(tmp$faqu_fix[cluster_i], na.rm=T)
    tmp$pol_fix[cluster_i] <- mean(tmp$fpol_fix[cluster_i], na.rm=T)
  }

  # independent SNPs (not in any LD cluster): own values, own cluster number
  indep_i <- which(is.na(tmp$cluster))
  if(length(indep_i) > 0){
    last_cluster <- max(c(0, tmp$cluster), na.rm=T)
    tmp$aq_fix[indep_i] <- tmp$faqu_fix[indep_i]
    tmp$pol_fix[indep_i] <- tmp$fpol_fix[indep_i]
    tmp$cluster[indep_i] <- (last_cluster + 1):(last_cluster + length(indep_i))
  }

  # a cluster whose markers are all NA would give NaN: fail rather than write it
  stopifnot(!anyNA(tmp$aq_fix), !anyNA(tmp$pol_fix))

  ### change fixation indexes to probability of drawing an m10 mutation
  tmp$pm10_aq <- tmp$aq_fix
  tmp$pm10_pol <- 1 - tmp$pol_fix

  # save per-chromosome AIMs file
  file_tag <- paste("AIMs_ch", chr, ".txt", sep="")

  write.table(tmp[, c("chr_id", "pos_0based", "cluster", "pm10_aq", "pm10_pol")],
              file_tag,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  print(paste("Finished file ", chr, sep=""))
}
