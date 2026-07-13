library(data.table)
library(igraph)
library(parallel)
library(LDscnR)

# ----------------------------
# LD edge precomputation
# ----------------------------

precompute_LD_edges <- function(
    GTs,
    map,
    r2_min = 0.1,
    max_bp = Inf,
    cores = 1
) {
  markers <- intersect(colnames(GTs), map$marker)
  if (length(markers) == 0) stop("No overlapping markers between GTs and map.")

  map_sub <- copy(map[marker %in% markers])
  setkey(map_sub, Chr, Pos)

  chr_levels <- unique(map_sub$Chr)

  out <- mclapply(chr_levels, function(ch) {
    chr_map <- map_sub[Chr == ch]
    chr_markers <- chr_map$marker

    if (length(chr_markers) < 2) {
      return(data.table(
        Chr = ch,
        marker1 = chr_markers,
        marker2 = chr_markers,
        r2 = 1,
        dist_bp = 0
      ))
    }

    gts <- as.matrix(GTs[, chr_markers, drop = FALSE])
    storage.mode(gts) <- "double"

    R2 <- cor(gts, use = "pairwise.complete.obs")^2
    R2[is.na(R2)] <- 0
    diag(R2) <- 0

    idx <- which(R2 >= r2_min, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

    if (nrow(idx) == 0) {
      return(data.table(
        Chr = ch,
        marker1 = chr_markers,
        marker2 = chr_markers,
        r2 = 1,
        dist_bp = 0
      ))
    }

    pos <- chr_map$Pos

    dt <- data.table(
      Chr = ch,
      marker1 = chr_markers[idx[, 1]],
      marker2 = chr_markers[idx[, 2]],
      r2 = R2[idx],
      dist_bp = abs(pos[idx[, 1]] - pos[idx[, 2]])
    )

    if (is.finite(max_bp)) dt <- dt[dist_bp <= max_bp]

    dt
  }, mc.cores = cores)

  out <- rbindlist(out, use.names = TRUE, fill = TRUE)
  setkey(out, Chr, marker1, marker2)

  out
}


# ----------------------------
# LD clusters from edge list
# ----------------------------
LD_igraph_components <- function(
    el,
    markers,
    r2_th = 0.8,
    bp_th = Inf
) {
  markers <- unique(markers)

  if (length(markers) == 0) {
    return(data.table(marker = character(), CL_id = integer(), n_loci = integer()))
  }

  if (length(markers) == 1) {
    return(data.table(marker = markers, CL_id = 1L, n_loci = 1L))
  }

  edges <- el[
    marker1 %in% markers &
      marker2 %in% markers &
      r2 >= r2_th
  ]

  if (is.finite(bp_th)) edges <- edges[dist_bp <= bp_th]

  if (nrow(edges) == 0) {
    return(data.table(
      marker = markers,
      CL_id = seq_along(markers),
      n_loci = 1L
    ))
  }

  g <- graph_from_data_frame(
    edges[, .(from = marker1, to = marker2)],
    directed = FALSE,
    vertices = data.table(name = markers)
  )

  comp <- components(g)

  clusters <- data.table(
    marker = names(comp$membership),
    CL_id = as.integer(comp$membership)
  )

  clusters[, n_loci := comp$csize[CL_id]]
  clusters
}


# ----------------------------
# Empty result, generalized
# ----------------------------
empty_result <- function(
    rho,
    th_ldw,
    th_XtX,
    p_names,
    r2_grid,
    lmin_grid
) {
  out <- CJ(r2_th = r2_grid, l_min = lmin_grid)

  out[, `:=`(
    th_ldw = th_ldw,
    th_XtX = th_XtX,
    rho = rho
  )]

  for (nm in p_names) {
    out[, (nm) := list(list(character()))]
  }

  out[]
}


# ----------------------------
# One grid point, generalized
# ----------------------------
run_one_grid <- function(
    map,
    el=NULL,
    ld_ws,
    rho,
    th_ldw,
    th_XtX,
    p_xtx_col,
    p_cols,
    p_names = names(p_cols),
    alpha = 0.05,
    r2_grid,
    lmin_grid,
    bp_th = Inf,
    cores
) {
  stopifnot(length(p_cols) == length(p_names))

  # Important: align ld_ws and map by marker before filtering
  common_markers <- intersect(map$marker, rownames(ld_ws))


  map_sub <- copy(map[marker %in% common_markers])
  ld_sub <- ld_ws[map_sub$marker, , drop = FALSE]

  keep_ld_w <- ld_sub[, rho] > quantile(ld_sub[, rho], th_ldw, na.rm = TRUE)


  #th_XtX = 0.9
  keep_XtX <- map_sub[[p_xtx_col]] < quantile(map_sub[[p_xtx_col]], 1 - th_XtX, na.rm = TRUE)
  #table(keep_ld_w)

  keep <- keep_ld_w & keep_XtX
  keep[is.na(keep)] <- FALSE
  #table(keep)

  if (!any(keep)) {
    return(cbind(empty_result(rho, th_ldw, th_XtX, p_names, r2_grid, lmin_grid),n_loci=length(which(keep))))
  }

  markers_keep <- map_sub[keep, marker]

  outliers <- setNames(vector("list", length(p_cols)), p_names)


  #i <- 1
  for (i in seq_along(p_cols)) {
    p_col <- p_cols[i]
    nm <- p_names[i]

    q <- p.adjust(unlist(map_sub[keep, ..p_col]), method = "fdr")
    outliers[[nm]] <- markers_keep[q < alpha]
  }

  if (length(unique(unlist(outliers))) == 0) {
    return(cbind(empty_result(rho, th_ldw, th_XtX, p_names, r2_grid, lmin_grid),n_loci=length(which(keep))))
  }


  if(is.null(el)){
    all_outliers <- unique(unlist(outliers))
    el <- precompute_LD_edges(
      GTs = GTs[, all_outliers, drop = FALSE],
      map = map_sub[marker %in% all_outliers],
      r2_min = 0.1,
      max_bp = 1e6,
      cores = 1
    )
  }

  out <- rbindlist(mclapply(r2_grid, function(r2_th) {
    clusters <- lapply(outliers, function(markers) {
      LD_igraph_components(
        el = el,
        markers = markers,
        r2_th = r2_th,
        bp_th = bp_th
      )
    })

    rbindlist(lapply(lmin_grid, function(l_min) {
      row <- data.table(
        r2_th = r2_th,
        l_min = l_min,
        th_ldw = th_ldw,
        th_XtX = th_XtX,
        rho = rho
      )

      for (nm in p_names) {
        row[, (nm) := list(clusters[[nm]][n_loci >= l_min, marker])]
      }

      row
    }), fill = TRUE)
  },mc.cores=cores), fill = TRUE)
  out[,n_loci:=length(which(keep))]

}


# ----------------------------
# Potential outliers, generalized
# ----------------------------
get_potential_outliers <- function(
    map,
    ld_ws,
    qt_grid,
    p_xtx_col,
    p_cols,
    alpha = 0.05
) {
  common_markers <- intersect(map$marker, rownames(ld_ws))

  #map <- copy(map[marker %in% common_markers])
  #ld <- ld_ws[map$marker, , drop = FALSE]

  potential <- character()

  for (rho in colnames(ld_ws)) {
    message("processing rho = ",rho)
    ld_vec <- ld_ws[, rho]

    potential <- unique(potential)
    for (th_ldw in qt_grid) {
      keep_ld_w <- ld_vec > quantile(ld_vec, th_ldw, na.rm = TRUE)

      for (th_XtX in qt_grid) {
        keep_XtX <- map[[p_xtx_col]] < quantile(map[[p_xtx_col]], 1 - th_XtX, na.rm = TRUE)

        keep <- keep_ld_w & keep_XtX
        keep[is.na(keep)] <- FALSE

        if (!any(keep)) next
        ##p_col <- p_cols[1]
        for (p_col in p_cols) {
          q <- p.adjust(unlist(map[keep, ..p_col]), method = "fdr")
          potential <- c(potential, map[keep, marker][q < alpha])
        }
      }
    }
  }
  return(unique(potential))

}

# ----------------------------
# Summerize
# ----------------------------
summarise_stability <- function(outliers, map, p_names) {
  map_C <- copy(map)

  for (nm in p_names) {
    C <- outliers[, table(unlist(get(nm))) / .N]
    C <- data.table(C)
    setnames(C, c("V1", "N"), c("marker", paste0("C_", nm)))

    map_C <- C[map_C, on = "marker"]
    map_C[is.na(get(paste0("C_", nm))), (paste0("C_", nm)) := 0]
  }

  map_C[]
}

# Example use
# ----------------------------
# read in raw data
# ----------------------------

if(!file.exists("data_1mb.rds")){


  data <- readRDS("./data/filtered_DIEM_1mb.rds")
  map_1mb <- data$map
  GTs_1mb <- data$GTs
  rm(data)
  gc()
  colnames(GTs_1mb) <- map_1mb$marker

  data <- readRDS("./data/filtered_DIEM_1mb_with_parents_MAF001.rds")
  sample_info <- data$sample_info

  ## read and transform BF to a pseudo-F scale
  ## from ./R/baypass.R
  BP_PC1 <- fread("./out_baypass/PC1_summary_betai_reg.out")
  map_1mb[, p_BF_PC1 := BP_PC1[,1 - rank(`BF(dB)`, ties.method = "average") / (.N + 1)]]

  BP_PC2 <- fread("./out_baypass/PC2_summary_betai_reg.out")
  map_1mb[, p_BF_PC2 := BP_PC2[,1 - rank(`BF(dB)`, ties.method = "average") / (.N + 1)]]

  xtx_PC1 <- fread("./out_baypass/PC1_DIEM_summary_pi_xtx.out")$XtXst
  map_1mb[, p_xtx_PC1 := 1 - rank(xtx_PC1, ties.method = "average") / (.N + 1)]


  library(GGally)
  dt <- data.table(
    XtX = xtx_PC1,
    ld_w_99 = ld_ws_1mb[, "0.99"],
    PC1 = BP_PC1$`BF(dB)`,
    PC2 = BP_PC2$`BF(dB)`
  )

  p1 <- ggpairs(
    dt,
    lower = list(
      continuous = wrap(
        "points",
        alpha = 0.3,
        size = 0.3
      )
    ),
    upper = list(
      continuous = wrap(
        "cor",
        size = 4
      )
    ),
    diag = list(
      continuous = wrap(
        "densityDiag",
        alpha = 0.5
      )
    )
  ) +
    theme_bw()

  png("./pairs.png",width = 8,height = 8,res = 100,units = "in")
  p1
  dev.off()
  map_1mb[, p_xtx_PC1 := 1 - rank(xtx_PC1, ties.method = "average") / (.N + 1)]

  gds <- create_gds_from_geno(geno=GTs_1mb[sample_info[,which(Species=="hybrid" & Population != "Sielva")],], map_1mb, "gds_formicia")

  #snpgdsClose(gds)
  ld_decay_corr <- compute_LD_decay(
    gds,
    el_data_folder = "./EL_corr_1mb/",
    ## for LD-decay and bg
    q = 0.95,
    ## for bg
    n_sub_bg = 5000,
    ## for decay
    n_win_decay = 20,
    max_pairs = 5000,
    max_SNPs_decay = Inf,
    n_strata = 20,
    overlap = 0.5,
    prob_robust = 0.95,
    keep_el = TRUE,
    slide=400, ## slide~400 is needed
    cores = 10,ld_method = "corr"
  )


  ld_ws_1mb <- precalculate_ld_w(rho=c(seq(0.25,0.95,by=0.05),0.99),ld_decay_corr)
  rownames(ld_ws_1mb) <- map_1mb$marker
  colnames(GTs_1mb) <- map_1mb$marker

  saveRDS(list(ld_decay=ld_decay_corr,ld_ws=ld_ws_1mb,sample_info=sample_info,map=map_1mb,GTs=GTs_1mb),"data_1mb.rds")
}else{
  data <- readRDS("data_1mb.rds")
  ld_decay_corr <- data$ld_decay
  data$ld_decay <- NULL
  ld_ws_1mb <- data$ld_ws
  data$ld_ws <- NULL
  sample_info <- data$sample_info
  data$sample_info <- NULL
  map_1mb <- data$map
  data$map <- NULL
  GTs_1mb <- data$GTs
  data$GTs <- NULL
  rm(data)
}
# analyse data ----------
#######################


r2_grid   <- seq(0.3, 0.8, by = 0.05)
lmin_grid <- (1:5)^2
qt_grid   <- c(seq(0.6, 0.95, by = 0.05),0.99)


p_cols <- c(
  PC1 = "p_BF_PC1",
  PC2 = "p_BF_PC2"
)

potential_outliers <- get_potential_outliers(
  map = map_1mb,
  ld_ws = ld_ws_1mb,
  qt_grid = qt_grid,
  p_xtx_col = "p_xtx_PC1",
  p_cols = p_cols,
  alpha = 0.05
)

GTs_pot <- GTs_1mb[, potential_outliers, drop = FALSE]
map_pot <- map_1mb[marker %in% potential_outliers]
ld_pot  <- ld_ws_1mb[map_pot$marker, , drop = FALSE]

el_potential <- precompute_LD_edges(
  GTs = GTs_1mb[, potential_outliers, drop = FALSE],
  map = map_1mb[marker %in% potential_outliers],
  r2_min = 0.1,
  max_bp = 1e6,
  cores = 4
)

param_grid <- CJ(
  rho = colnames(ld_ws_1mb),
  th_ldw = qt_grid,
  th_XtX = qt_grid
)
#i <- 1
outliers <- rbindlist(
  mclapply(seq_len(nrow(param_grid)), function(i) {
    pars <- param_grid[i]
    cat(i,"..")
    out <- run_one_grid(
      map = map_1mb,
      el = el_potential,
      ld_ws = ld_ws_1mb,
      rho = pars$rho,
      th_ldw = pars$th_ldw,
      th_XtX = pars$th_XtX,
      p_xtx_col = "p_xtx_PC1",
      p_cols = p_cols,
      alpha = 0.05,
      r2_grid = r2_grid,
      lmin_grid = lmin_grid,
      bp_th = Inf,
      cores = 1
    )
  },mc.cores=1),
  fill = TRUE
)
#save.image()
#gc()
#saveRDS(outliers,"outliers_1200.rds")
outliers[,plot(as.numeric(rho),lengths(PC2),col=as.numeric(factor(th_XtX)))]


outliers[th_ldw>0.6 & th_XtX>0.6]

map_C <- summarise_stability(
  outliers = outliers[th_XtX>=0.8 & l_min>2 & rho>=0.6],
  map = map_1mb,
  p_names = names(p_cols)
)

par(mfcol=c(2,1))
map_C[,plot(C_PC2)]
map_C[,plot(C_PC1)]

joint_outliers <- map_C[C_PC2>0.05 | C_PC1>0.05,marker]

el_joint_outl <- precompute_LD_edges(
  GTs = GTs_1mb[, joint_outliers, drop = FALSE],
  map = map_1mb[marker %in% joint_outliers],
  r2_min = 0.1,
  max_bp = 1e6,
  cores = cores
)

final_ORs <- LD_igraph_components(bp_th = 1e6,r2_th = 0.25,markers = joint_outliers,el = el_joint_outl)

OR_cls <- split(final_ORs$marker,final_ORs$CL_id)
