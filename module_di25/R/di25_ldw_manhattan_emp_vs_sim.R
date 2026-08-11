## =============================================================================
## module_di25 -- local-LD (ld_w_095) windowed-median Manhattan: empirical vs SIM
## =============================================================================
## Per-SNP local LD ld_w_095 = median composite r^2 among neighbours within the
## rho=0.95 LD-decay distance (compute_ld_w), computed with the SAME DI25 pipeline
## on both sides (LD-decay hybrids-only, slide=100, ld_method="corr"):
##   * empirical -- cached DI25 ld_w_095  (module_di25/data/di25_stage1.rds)
##   * simulated -- ONE DIEM bootstrap replicate (data/diem_outs_demo/diem_boot<N>_output.bed)
##
## Both use the same ~51.6k DI25 diagnostic markers and the 165 hybrids only
## (parents dropped -- they are near-fixed at diagnostic markers and manufacture
## artificial LD). r^2 is orientation-invariant, so the DIEM 0/1/2 states are used
## directly as dosage without re-polarising.
##
## The figure is the windowed MEDIAN of ld_w_095 along the genome, faceted by
## chromosome, empirical vs simulated overlaid -- the ld_w analogue of
## Figures_legacy/moduleE_ldw_manhattan_foundersweep.pdf.
##
## Run from the formica_hybrid repo root:
##   Rscript module_di25/R/di25_ldw_manhattan_emp_vs_sim.R [rep]   (rep default 1)
## =============================================================================
suppressMessages({
  library(data.table); library(ggplot2); library(SNPRelate)
})
devtools::load_all("~/gitlab/LDscnR/")

args   <- commandArgs(trailingOnly = TRUE)
REP    <- if (length(args) >= 1) as.integer(args[1]) else 1L
STAGE1 <- "module_di25/data/di25_stage1.rds"
BED    <- sprintf("data/diem_outs_demo/diem_boot%d_output.bed", REP)
GDS    <- sprintf("module_di25/data/di25_sim_rep%d.gds", REP)
OUTRDS <- sprintf("module_di25/data/di25_ldw_emp_vs_sim_rep%d.rds", REP)
OUTPNG <- sprintf("module_di25/Figures/di25_ldw_manhattan_emp_vs_sim_rep%d.png", REP)
OUTPDF <- sub("\\.png$", ".pdf", OUTPNG)

RHO_LDW <- 0.95
CORES   <- 4
WIN_BP  <- 500000L   # windowed-median bin size
MIN_N   <- 5L        # min markers per window to report a median

## ---- empirical DI25 ld_w_095 (cached) --------------------------------------
message("[emp] loading cached DI25 ld_w_095  <- ", STAGE1)
emp <- as.data.table(readRDS(STAGE1)$map)[, .(Chr, Pos, ld_w = ld_w_095)]
emp[, chr_id := as.integer(sub("Chr", "", Chr))]

## ---- simulated: parse one DIEM bed, hybrids only, DI25 pipeline ld_w --------
message("[sim] parsing ", BED)
hdr  <- readLines(BED, n = 2)
inds <- strsplit(strsplit(hdr[2], "\t")[[1]][10], "\\|")[[1]]
sim  <- fread(BED, skip = 2, header = FALSE, sep = "\t",
              select = c(1, 3, 10), colClasses = list(character = c(1, 10)))
chr_id <- as.integer(sub("ch", "", sim$V1))
pos    <- as.integer(sim$V3)
geno   <- sub("^S", "", sim$V10)

## explode DIEM state strings -> markers x individuals dosage (0/1/2, '_'/'.'=NA)
S <- matrix(unlist(strsplit(geno, "", fixed = TRUE), use.names = FALSE),
            nrow = length(geno), byrow = TRUE)
stopifnot(ncol(S) == length(inds))
dos <- matrix(NA_integer_, nrow(S), ncol(S))
dos[S == "0"] <- 0L; dos[S == "1"] <- 1L; dos[S == "2"] <- 2L
colnames(dos) <- inds

## hybrids only for LD (drop aq_/pol_ parents), map ordered by (chr, pos)
hyb  <- grep("^hyb_", inds)
map  <- data.table(Chr = paste0("Chr", chr_id), Pos = pos,
                   marker = paste0("Chr", chr_id, ":", pos), chr_id = chr_id)
ord  <- order(map$chr_id, map$Pos)
map  <- map[ord]; dos <- dos[ord, , drop = FALSE]
GTs_hyb <- t(dos[, hyb, drop = FALSE])                     # hybrids x markers
message(sprintf("  %d markers | %d hybrids (LD)", nrow(map), nrow(GTs_hyb)))

message("[sim] LD decay (hybrids only) + ld_w rho=", RHO_LDW)
gds <- create_gds_from_geno(geno = GTs_hyb, map = map, gds_path = GDS)
ld_decay <- compute_LD_decay(gds, keep_el = TRUE, slide = 100,
                             ld_method = "corr", cores = CORES)
snpgdsClose(gds)
lw  <- compute_ld_w(ld_decay, rho = RHO_LDW, cores = CORES)
ids <- unlist(lapply(ld_decay$by_chr, function(o) o$snp_ids), use.names = FALSE)
map[, ld_w := lw[match(marker, ids)]]
sim_dt <- map[, .(Chr, Pos, ld_w, chr_id)]

## ---- windowed median per (chromosome, source) ------------------------------
win_median <- function(dt, src) {
  d <- dt[is.finite(ld_w)]
  d[, win := (Pos %/% WIN_BP)]
  d[, .(mid = (unique(win) * WIN_BP + WIN_BP / 2), n = .N,
        ld_w_med = median(ld_w, na.rm = TRUE)),
    by = .(chr_id, Chr, win)][n >= MIN_N][, src := src][]
}
wm <- rbind(win_median(emp, "empirical"), win_median(sim_dt, "simulated (rep)"))
wm[, Chr := factor(Chr, levels = paste0("Chr", sort(unique(wm$chr_id))))]
wm[, src := factor(src, levels = c("empirical", "simulated (rep)"))]

saveRDS(list(wm = wm, emp = emp, sim = sim_dt, rep = REP,
             win_bp = WIN_BP, rho = RHO_LDW), OUTRDS)

## ---- plot: windowed-median ld_w Manhattan, faceted by chromosome -----------
cols <- c("empirical" = "#d95f02", "simulated (rep)" = "#1b9e77")
p <- ggplot(wm, aes(mid / 1e6, ld_w_med, colour = src)) +
  geom_line(linewidth = 0.4) +
  facet_grid(. ~ Chr, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "position (Mbp)", y = expression("windowed median local LD  " * ld[w][0.95]),
       title = sprintf("DI25 local-LD (ld_w_0.95) landscape: empirical vs simulated replicate %d", REP),
       subtitle = sprintf("median composite r^2 (rho=%.2f), hybrids only, %dkb windows (>=%d markers); same DI25 marker panel",
                          RHO_LDW, WIN_BP %/% 1000L, MIN_N)) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.spacing = unit(1, "pt"),
        strip.text.x = element_text(size = 5, angle = 90),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top")
ggsave(OUTPDF, p, width = 14, height = 4.2)
ggsave(OUTPNG, p, width = 14, height = 4.2, dpi = 150)

## ---- summary: level + per-chromosome landscape Spearman --------------------
lvl <- wm[, .(mean_win_med = mean(ld_w_med, na.rm = TRUE)), by = src]
pc  <- dcast(wm, chr_id ~ src, value.var = "ld_w_med", fun.aggregate = mean)
rho_land <- cor(pc[["empirical"]], pc[["simulated (rep)"]], method = "spearman", use = "complete.obs")
cat("\n=== ld_w_0.95 windowed-median summary ===\n")
print(lvl)
cat(sprintf("per-chromosome landscape Spearman (emp vs sim): %.3f\n", rho_land))
cat("saved: ", OUTRDS, "\n       ", OUTPNG, "\n       ", OUTPDF, "\n")
