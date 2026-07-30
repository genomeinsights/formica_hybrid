## =============================================================================
## Compare EMPIRICAL hybrids vs COLLEAGUE'S SIMULATED hybrids (DIEM-parsed)
## =============================================================================
## A collaborator produced a simulated replicate of the 20-population hybrid
## design and parsed it through DIEM exactly as the empirical data were parsed
## (R/LD_decay_sim.R -> data/sim_parsed_data/tst.rds). Her SNP set is slightly
## different from ours, so EVERYTHING here is restricted to the INTERSECTION of
## markers between the two datasets.
##
## Both datasets have the same structure: 20 hybrid demes (sizes 3-20) + aq + pol
## parents. We compare, on the shared markers, the two properties the biology
## turns on:
##   (1) SORTING   -- directional ancestry fixation (parallelism_stats/sort_class),
##                    stratified by (empirical) DiagnosticIndex -> dose-response.
##   (2) LD        -- pooled composite-r^2 decay with physical distance
##                    (compute_LD_decay, colleague's exact conventions).
##
## Marker-matching note: the sim map labels chromosomes "ch1", the empirical set
## "Chr1"; positions are on the same reference. We normalise ch->Chr and match on
## "Chr:Pos". Result: 51,065 shared markers (95.8% of the sim set). The shared set
## is all moderate/diagnostic DI (emp DI -38..-4); it contains NO strictly-neutral
## markers, so there is no internal neutral control here -- this is a like-for-like
## empirical-vs-sim comparison, not a null test.
##
## Design decisions (see inline):
##  - Sorting is GATED and ORIENTED by the EMPIRICAL parents for BOTH datasets, so
##    the marker set, the "aquilonia allele", and the maf gate are identical on
##    both sides; only the hybrid dosages differ. (A sim-own-parent orientation is
##    also computed as a robustness check.)
##  - LD is POOLED across all demes for both (Wahlund/among-pop LD is the intended
##    diagnostic, per the locked convention), same compute_LD_decay settings.
##
## STAGE env var: "sort" (fast), "ld" (heavy), or "all" (default).
## Outputs: data/moduleE_sim/compare_emp_vs_simDIEM.rds
##          Figures/compare_emp_vs_simDIEM_sorting.pdf / _ld.pdf
## =============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork); library(LDscnR) })

FORMICA <- "/Users/petrikem/gitlab/formica_hybrid"
SIM     <- file.path(FORMICA, "data/sim_parsed_data/tst.rds")
BUNDLE  <- file.path(FORMICA, "moduleE_slim/inputs/empirical_bundle.rds")
EMP_RD  <- file.path(FORMICA, "data/hybrids_only_maf005.Rdata")
OUTDIR  <- file.path(FORMICA, "data/moduleE_sim"); dir.create(OUTDIR, showWarnings=FALSE)
FIGDIR  <- file.path(FORMICA, "Figures");         dir.create(FIGDIR, showWarnings=FALSE)
GDSDIR  <- file.path(OUTDIR, "gds_compare");       dir.create(GDSDIR, showWarnings=FALSE)
source(file.path(FORMICA, "moduleE_slim/scripts/parallelism_stats.R"))
STAGE <- Sys.getenv("STAGE", "all")

## ---- shared marker set ------------------------------------------------------
sim  <- readRDS(SIM)
smap <- as.data.table(sim$map)
smap[, mk := paste0(sub("^ch","Chr",Chr), ":", Pos)]      # normalise ch1 -> Chr1
b    <- readRDS(BUNDLE)
uni  <- as.data.table(b$universe); setkey(uni, marker)
shared <- intersect(smap$mk, uni$marker)
shared <- uni[shared][order(Chr, Pos), marker]            # order by genome position
cat(sprintf("shared markers: %d (%.1f%% of sim's %d)\n",
            length(shared), 100*length(shared)/nrow(smap), nrow(smap)))

DIemp <- setNames(uni[shared, DI],         shared)        # empirical DI (stratifier)
pmaf  <- setNames(uni[shared, parent_maf], shared)        # empirical pooled-parent maf (gate)
f_aq  <- b$f_aq_par[shared]; f_pol <- b$f_pol_par[shared] # empirical parental freqs (orient)

## DI bins (locked Module A / HANDOFF bins, restricted to the range present here)
BINS <- c(-Inf,-40,-30,-25,-20,-15,-10,Inf)
dibin <- cut(DIemp, BINS)

## =============================================================================
## (1) SORTING
## =============================================================================
## helper: per-marker sort_class for a pops x markers DOSAGE matrix, oriented &
## gated by the empirical parents (common frame for both datasets).
run_sort <- function(mm, hyb_pops){
  pm <- rbind(mm, aquilonia_parent = f_aq*2, polyctena_parent = f_pol*2)
  colnames(pm) <- shared
  r <- parallelism_stats(list(pop_means = pm),
                         hybrid_pops = hyb_pops,
                         aqu_pops = "aquilonia_parent", pol_pops = "polyctena_parent",
                         orient = "parents",
                         parent_maf = pmaf, min_parent_maf = 0.15,
                         fix_th = 0.15, sort_th = 0.5, null_prob = 0.5,
                         DI = DIemp)
  r[match(shared, marker)]                                # returns a data.table keyed by marker
}

## --- EMPIRICAL: per-pop dosage from the compact bundle -----------------------
emp_M <- b$emp_mean[, shared, drop=FALSE] / 1000          # pops x markers dosage (0-2)
sort_emp <- run_sort(emp_M, rownames(emp_M))

## --- SIM: per-pop dosage from the parsed genotypes ---------------------------
GT <- sim$GTs
colnames(GT) <- smap$mk[match(colnames(GT), smap$marker)] # -> Chr:Pos names
GT <- GT[, shared, drop=FALSE]
is_aq  <- grepl("^aq",  rownames(GT)); is_pol <- grepl("^pol", rownames(GT))
is_hyb <- !(is_aq | is_pol)
hybpop <- sub("_.*", "", rownames(GT)[is_hyb])            # hyb11, hyb115, ...
sim_pops <- sort(unique(hybpop))
sim_M <- t(vapply(sim_pops, function(p)
  colMeans(GT[is_hyb,,drop=FALSE][hybpop==p,,drop=FALSE], na.rm=TRUE),
  numeric(length(shared))))
sort_sim <- run_sort(sim_M, sim_pops)

## robustness: SIM oriented by its OWN parents (guards against coding flips) ---
f_aq_sim  <- colMeans(GT[is_aq,,drop=FALSE],  na.rm=TRUE)/2
f_pol_sim <- colMeans(GT[is_pol,,drop=FALSE], na.rm=TRUE)/2
run_sort_own <- function(mm, hyb_pops){
  pm <- rbind(mm, aquilonia_parent=f_aq_sim*2, polyctena_parent=f_pol_sim*2)
  colnames(pm) <- shared
  r <- parallelism_stats(list(pop_means=pm), hybrid_pops=hyb_pops,
        aqu_pops="aquilonia_parent", pol_pops="polyctena_parent", orient="parents",
        parent_maf=pmaf, min_parent_maf=0.15, fix_th=0.15, sort_th=0.5,
        null_prob=0.5, DI=DIemp)
  r[match(shared, marker)]
}
sort_sim_own <- run_sort_own(sim_M, sim_pops)

## --- Hudson pairwise F_ST helper (mean over pop pairs), per marker set -------
hudson_mean <- function(M, ns){                           # M: pops x markers dosage
  P <- M/2; K <- nrow(P); n2 <- ns*2
  num <- den <- 0
  for(i in 1:(K-1)) for(j in (i+1):K){
    p1<-P[i,];p2<-P[j,]
    nu<-(p1-p2)^2 - p1*(1-p1)/(n2[i]-1) - p2*(1-p2)/(n2[j]-1)
    de<-p1*(1-p2)+p2*(1-p1); ok<-is.finite(nu)&is.finite(de)&de>0
    num<-num+sum(nu[ok]); den<-den+sum(de[ok]) }
  max(num/den, 0)
}
emp_ns <- b$emp_ns[rownames(emp_M)]
sim_ns <- as.integer(table(hybpop)[sim_pops])

## --- summarise sorting by DI bin ---------------------------------------------
sorted_frac <- function(sc) mean(sc %in% c("aquilonia","polyctena"), na.rm=TRUE)
aq_lead     <- function(sc){ a<-sum(sc=="aquilonia",na.rm=TRUE); p<-sum(sc=="polyctena",na.rm=TRUE); (a-p)/max(a+p,1) }
sortsum <- rbindlist(lapply(levels(dibin), function(lv){
  ix <- which(dibin==lv)
  data.table(DI_bin=lv, n_markers=length(ix),
    emp_sorted = sorted_frac(sort_emp$sort_class[ix]),
    sim_sorted = sorted_frac(sort_sim$sort_class[ix]),
    emp_aqlead = aq_lead(sort_emp$sort_class[ix]),
    sim_aqlead = aq_lead(sort_sim$sort_class[ix]),
    emp_fst    = hudson_mean(emp_M[,ix,drop=FALSE], emp_ns),
    sim_fst    = hudson_mean(sim_M[,ix,drop=FALSE], sim_ns))
}))
cat("\n=== SORTING by empirical-DI bin (empirical-parent orientation, shared markers) ===\n")
print(sortsum)
cat(sprintf("\nsim orientation robustness: cor(emp-oriented sim sort dir, own-oriented) over classified markers:\n"))
both <- !is.na(sort_sim$sort_class) & !is.na(sort_sim_own$sort_class)
cat(sprintf("  agreement of sort_class = %.3f (n=%d)\n",
            mean(sort_sim$sort_class[both]==sort_sim_own$sort_class[both]), sum(both)))

## =============================================================================
## (2) LD DECAY (pooled across demes; colleague's compute_LD_decay settings)
## =============================================================================
LD_ARGS <- list(q=0.95, n_sub_bg=5000, n_win_decay=20, max_pairs=5000,
                max_SNPs_decay=Inf, n_strata=20, overlap=0.5, prob_robust=0.95,
                keep_el=TRUE, slide=1000, cores=8, ld_method="corr")
ld_of <- function(geno, tag){                             # geno: inds x shared markers
  mp <- uni[shared, .(Chr, Pos, marker)]
  gpath <- file.path(GDSDIR, paste0(tag, ".gds"))
  if (file.exists(gpath)) unlink(gpath)
  gds <- create_gds_from_geno(geno = geno, map = as.data.frame(mp), gds_path = gpath)
  do.call(compute_LD_decay, c(list(gds), LD_ARGS))
}

ld_emp <- ld_sim <- NULL
if (STAGE %in% c("ld","all")){
  cat("\n=== LD decay on shared markers (this is the heavy step) ===\n")
  e <- new.env(); load(EMP_RD, envir=e)
  geno_emp <- e$GTs_hybrids_005[, shared, drop=FALSE]      # 164 hybrids x shared
  rm(e); gc()
  ld_emp <- ld_of(geno_emp, "emp_shared")
  geno_sim <- GT[is_hyb, , drop=FALSE]                     # 165 hybrids x shared
  ld_sim <- ld_of(geno_sim, "sim_shared")
  ## scalar decay summary: median bp to decay to r=0.95 (per-chr), pooled model a
  ds <- function(ld) c(median_slide_bp_095 = median(ld$decay_sum$slide_bp_rho_0.95, na.rm=TRUE),
                       median_slide_snp_095= median(ld$decay_sum$slide_snp_rho_0.95, na.rm=TRUE))
  cat("EMP decay:", paste(names(ds(ld_emp)), round(ds(ld_emp)), collapse="  "), "\n")
  cat("SIM decay:", paste(names(ds(ld_sim)), round(ds(ld_sim)), collapse="  "), "\n")
}

## =============================================================================
## save + figures
## =============================================================================
saveRDS(list(shared=shared, DI=DIemp, sortsum=sortsum,
             sort_emp=sort_emp, sort_sim=sort_sim, sort_sim_own=sort_sim_own,
             ld_emp=ld_emp, ld_sim=ld_sim,
             emp_pops=rownames(emp_M), sim_pops=sim_pops,
             emp_ns=emp_ns, sim_ns=sim_ns),
        file.path(OUTDIR, "compare_emp_vs_simDIEM.rds"))

## sorting dose-response figure
sl <- melt(sortsum, id.vars=c("DI_bin","n_markers"),
           measure.vars=patterns("_sorted$|_fst$|_aqlead$"))
sl[, c("src","stat") := tstrsplit(variable, "_", fixed=TRUE)]
sl[, DI_bin := factor(DI_bin, levels=levels(dibin))]
labmap <- c(sorted="sorted fraction", fst="Hudson F_ST", aqlead="aq-lead (dir. bias)")
p_sort <- ggplot(sl, aes(DI_bin, value, colour=src, group=src)) +
  geom_line() + geom_point(size=1.6) +
  facet_wrap(~factor(labmap[stat], levels=labmap), scales="free_y") +
  scale_colour_manual(values=c(emp="#d95f02", sim="#1b9e77"), name=NULL,
                      labels=c(emp="empirical", sim="simulated")) +
  labs(x="empirical DiagnosticIndex bin", y=NULL,
       title="Empirical vs simulated hybrids: sorting & differentiation by DI (shared markers)") +
  theme_bw(base_size=8) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(FIGDIR,"compare_emp_vs_simDIEM_sorting.pdf"), p_sort, width=10, height=3.6)
ggsave(file.path(FIGDIR,"compare_emp_vs_simDIEM_sorting.png"), p_sort, width=10, height=3.6, dpi=130)

## LD decay figure (raw binned r^2 vs distance from kept edge lists)
if (!is.null(ld_emp) && !is.null(ld_sim)){
  bin_el <- function(ld, tag){                            # edge-list cols: d (bp), r2
    el <- rbindlist(lapply(ld$by_chr, function(z) if(!is.null(z$el)) z$el[, .(d, r2)]), fill=TRUE)
    if(!nrow(el)) return(NULL)
    el <- el[is.finite(d) & is.finite(r2) & d > 0]
    el[, db := cut(d, breaks=c(0,10^seq(2,7,0.25)), labels=FALSE)]
    el[, .(dist=as.double(median(d)), ld=mean(r2), n=.N, src=tag), by=db][!is.na(db)][order(dist)]
  }
  be <- rbind(bin_el(ld_emp,"empirical"), bin_el(ld_sim,"simulated"))
  if(!is.null(be) && nrow(be)){
    p_ld <- ggplot(be, aes(dist, ld, colour=src)) + geom_line() + geom_point(size=1) +
      scale_x_log10() + scale_colour_manual(values=c(empirical="#d95f02", simulated="#1b9e77"), name=NULL) +
      labs(x="distance (bp, log)", y="mean composite r^2",
           title="Pooled LD decay: empirical vs simulated (shared markers)") +
      theme_bw(base_size=9)
    ggsave(file.path(FIGDIR,"compare_emp_vs_simDIEM_ld.pdf"), p_ld, width=6, height=4)
    ggsave(file.path(FIGDIR,"compare_emp_vs_simDIEM_ld.png"), p_ld, width=6, height=4, dpi=130)
  }
}
cat("\nsaved: data/moduleE_sim/compare_emp_vs_simDIEM.rds + Figures/compare_emp_vs_simDIEM_*\n")
