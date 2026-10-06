## Shared helpers for module_external_validation EMMAX scripts (source after devtools::load_all("LDscnR")).

## mean-impute missing dosages per SNP
meanimp <- function(M) { mu <- colMeans(M, na.rm = TRUE); w <- which(is.na(M), arr.ind = TRUE); M[w] <- mu[w[, 2]]; M }

## GLS effect, SE and variance components for selected columns of an emmax_setup() object; mirrors emmax_fast()
beta_se <- function(prep, y, idx) {
  re <- emma.REMLE(y, prep$Xo, prep$Kn, eig.R = prep$eigR)
  wv <- 1 / sqrt(re$vg * prep$lam + re$ve)
  yt <- as.numeric(crossprod(prep$V, y)) * wv; xo <- prep$xot * wv; Xtw <- prep$Xt[, idx, drop = FALSE] * wv
  a2 <- sum(xo * xo); ra <- yt - xo * (sum(xo * yt) / a2)
  Rb <- Xtw - outer(xo, as.numeric(crossprod(xo, Xtw)) / a2)
  b <- as.numeric(crossprod(ra, Rb)) / colSums(Rb * Rb)
  rss <- sum(ra * ra) - b^2 * colSums(Rb * Rb); se <- sqrt(rss / prep$df2 / colSums(Rb * Rb))
  list(beta = b, se = se, vg = re$vg, ve = re$ve, h2 = re$vg / (re$vg + re$ve))
}

## rank-based inverse-normal transform
rint <- function(y) qnorm((rank(y) - 0.5) / length(y))
