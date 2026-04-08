library(data.table)
#### --- parse and save GTs data --- ####
#devtools::install_github("genomeinsights/LDscnR")


## --- read in DIEM data --- ##
DIEM_2mb <- fread("./data/Formica_hybrid_populations_diem_output.bed")

## hybrid sample data
sample_info <- fread("./data/Sample_covariate_info_outlier_analysis_20.txt")
## only autosomal and bi-allelic
DIEM_2mb <- DIEM_2mb[`#Chrom`!="mtDNA" & nVNTs==2]


## sample names
DIEM_samples <- colnames(DIEM_2mb)[10]
DIEM_samples <- strsplit(DIEM_samples,"|",fixed = TRUE)[[1]]

## map
map_DIEM <- DIEM_2mb[,.(Chr=gsub("chromosome_","Chr",`#Chrom`),Pos=End,SNV,Qual,DiagnosticIndex,Support)]
map_DIEM[,marker := paste(Chr,Pos,sep=":")]

#new_order <- map_DIEM[,order(as.numeric(gsub("Chr","",Chr)),Pos)]

# Pre-extract only the column you need so workers do not carry full DIEM_2mb
track_strings <- as.character(DIEM_2mb[[10]])

# Process a chunk of rows per worker, not one row per future
idx_chunks <- split(seq_len(nrow(DIEM_2mb)),
                    ceiling(seq_len(nrow(DIEM_2mb)) / 1000))

pb <- txtProgressBar(min = 0, max = nrow(map_DIEM)-1, style = 3)

keep <- match(sample_info$Sample_ID,DIEM_samples,nomatch = 0)

n_inds <- length(keep)

DIEM_data <- rbindlist(lapply(idx_chunks,function(idx){
  setTxtProgressBar(pb, max(idx))
  out <- vector("list", length(idx))
  for (j in seq_along(idx)) {
    x <- idx[j]

    parts <- strsplit(track_strings[x], "|", fixed = FALSE)[[1]]
    track <- suppressWarnings(as.numeric(parts[-1]))[keep]

    nonmiss_keep <- sum(!is.na(track))

    if (nonmiss_keep == 0) {
      out[[j]] <- NULL
      next
    }

    #skip if too many missing
    missingness = sum(is.na(track)) / n_inds

    if (missingness >= 0.15) {
      out[[j]] <- NULL
      next
    }

    ## filter for maf
    maf <- sum(track, na.rm = TRUE) / nonmiss_keep / 2
    maf <- min(maf, 1 - maf)

    if (maf <= 0.1) {
      out[[j]] <- NULL
      next
    }

    out[[j]] <- cbind(map_DIEM[x], maf=maf,track=list(track))
  }
  rbindlist(out[!vapply(out, is.null, logical(1))])
}))

map_DIEM <- DIEM_data[,.(Chr,Pos,marker,SNV,Qual,DiagnosticIndex,Support)]
GTs_DIEM <- do.call(cbind,DIEM_data$track)
rownames(GTs_DIEM) <- sample_info$Sample_ID
colnames(GTs_DIEM) <- map_DIEM$marker
## save data ordered by chromosme and SNP
new_order <- map_DIEM[,order(as.numeric(gsub("Chr","",Chr)),Pos)]
saveRDS(list(GTs=GTs_DIEM[,new_order],map=map_DIEM[new_order],sample_info=sample_info),"./data/filtered_DIEM_1mb.rds")
q("no")
