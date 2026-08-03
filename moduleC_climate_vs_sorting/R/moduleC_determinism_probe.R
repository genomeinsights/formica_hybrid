## =========================================================
## MODULE C (revised) -- BayPass determinism probe
## =========================================================
## The mandatory faithful-regeneration gate (max|dk| == 0) assumes BayPass with
## `-nthreads 6 -seed 74` is bit-reproducible. If threaded reductions were
## non-deterministic, a few eMLGs whose null BF sits within floating-point epsilon
## of the observed BF could flip, giving a benign but non-zero dk that would only
## surface after the ~18 h full run. This probe tests determinism cheaply BEFORE
## the long run: it runs BayPass twice on the SAME small covariate file with the
## exact null options and asserts the two BF vectors are bit-identical.
##
## Reads : aland_excluded_eMLG/null/{u_eMLG.geno,omega_mat_omega.out,u_DIEM.size,null_b01.env}
## Writes: aland_excluded_eMLG/null/probe_{a,b}*  (deleted on success)
## Run from the formica_hybrid repo root. ~10-20 min.

suppressMessages(library(data.table))
ND      <- "aland_excluded_eMLG/null"
BAYPASS <- "/Users/petrikem/gitlab/baypass_public-master/sources/g_baypass"
OPT     <- "-nthreads 6 -nocovscaling -nval 500 -burnin 5000 -thin 25 -seed 74"
NPROBE  <- 3                                   # a few covariates suffice to test threading
stopifnot(file.exists(BAYPASS))

## small efile = first NPROBE null covariates from batch 1 (same for both runs)
Y  <- as.matrix(fread(file.path(ND, "null_b01.env")))   # 1000 rows (covariates) x 19 pops
ef <- file.path(ND, "probe.env")
write.table(Y[seq_len(NPROBE), , drop = FALSE], ef,
            quote = FALSE, row.names = FALSE, col.names = FALSE)

run_probe <- function(tag) {
  pref <- file.path(ND, paste0("probe_", tag))
  st <- system(paste0(BAYPASS, " -countdatafile ", ND, "/u_eMLG.geno -omegafile ", ND,
    "/omega_mat_omega.out -efile ", ef, " -poolsizefile ", ND, "/u_DIEM.size ",
    OPT, " -outprefix ", pref, " > ", pref, "_stdout.log 2>&1"))
  stopifnot("probe BayPass run failed" = st == 0)
  fread(paste0(pref, "_summary_betai_reg.out"))$`BF(dB)`
}

message("[probe] run A ...");  a <- run_probe("a")
message("[probe] run B ...");  b <- run_probe("b")
stopifnot("probe BF length mismatch" = length(a) == length(b))
d <- max(abs(a - b))
cat(sprintf("\n[probe] identical rows: %d/%d ; max|BF_A - BF_B| = %.3e\n",
            sum(a == b), length(a), d))
deterministic <- (d == 0)
if (deterministic) {
  cat("[probe] RESULT: BayPass is bit-deterministic under these options -> ",
      "the exact max|dk|==0 gate is safe.\n", sep = "")
} else {
  cat("[probe] RESULT: BayPass is NOT bit-deterministic (max diff ", d, "). ",
      "The strict max|dk|==0 gate may trip on benign boundary flips; revisit the ",
      "gate (e.g. allow dk within a documented small bound) before the full run.\n", sep = "")
}
## clean up
file.remove(Sys.glob(file.path(ND, "probe_*")))
file.remove(ef)
cat(sprintf("[probe] deterministic = %s\n", deterministic))
