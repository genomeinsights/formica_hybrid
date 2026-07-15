library(data.table)
library(LDscnR)

## LD-pruning/complexity-reduction for this simulated-data/EMMAX pipeline now
## lives in the package itself: LDscnR::ld_complexity_reduction(GTs, map,
## LD_decay, rho, cores, ld_w_col, ld_w_threshold). It also feeds CL_id/n_loci
## columns compatible with make_eMLGs(map_cl = ...).

files <- list.files("./data/diem_outs/",full.names = TRUE)
n_files <- 1
#file <- files[1]
lapply(files[seq_len(n_files)],function(file){
  boot_rep <- do.call(rbind,strsplit(file,"_"))[,4]
  
  DIEM <- fread(file)
  DIEM_samples <- colnames(DIEM)[10]
  DIEM_samples <- strsplit(DIEM_samples,"|",fixed = TRUE)[[1]]
  
  
  ## map
  map <- DIEM[,.(Chr=gsub("ch","Chr",`#Chrom`),Pos=End)]
  map[,marker := paste(Chr,Pos,sep=":")]
  
  # Pre-extract only the column you need so workers do not carry full DIEM
  track_strings <- as.character(DIEM[[10]])
  
  # Process a chunk of rows per worker, not one row per future
  idx_chunks <- split(seq_len(nrow(DIEM)),ceiling(seq_len(nrow(DIEM)) / 1000))
  
  pb <- txtProgressBar(min = 0, max = nrow(map)-1, style = 3)
  
  
  DIEM_data <- rbindlist(lapply(idx_chunks,function(idx){
    setTxtProgressBar(pb, max(idx))
    out <- vector("list", length(idx))
    for (j in seq_along(idx)) {
      x <- idx[j]
      
      parts <- strsplit(track_strings[x], "|", fixed = FALSE)[[1]]
      track <- suppressWarnings(as.numeric(parts[-1]))
      
      out[[j]] <- cbind(map[x], track=list(track))
    }
    rbindlist(out[!vapply(out, is.null, logical(1))])
  }))
  
  GTs <- do.call(cbind,DIEM_data$track)
  rownames(GTs) <- DIEM_samples
  colnames(GTs) <-  map$marker
  saveRDS(list(GTs=GTs,map=map),paste0("./data/sim_parsed_data/",boot_rep,".rds"))
  
})
