library(LDscnR)
library(ggplot2)
data <- readRDS("./data/sim_parsed_data/boot1.rds")
GTs_sim <- data$GTs
map_sim <- data$map
samples_sim <- rownames(GTs_sim)

hybrids_sim <- which(!grepl("pol",samples_sim) | grepl("aq",samples_sim))
gds_sim <- create_gds_from_geno(geno=GTs_sim[hybrids_sim,], map_sim, "sim.gds")
map_sim[,maf_hyb:= snpgdsSNPRateFreq(gds_sim)$MinorFreq]
#snpgdsClose(gds_sim)

ld_decay_sim <- compute_LD_decay(
  gds_sim,
  ## for LD-decay and bg
  q = 0.95,
  ## for bg
  n_sub_bg = 5000,
  ## for decay
  n_win_decay = 100,
  max_pairs = 5000,
  max_SNPs_decay = Inf,
  n_strata = 20,
  overlap = 0.5,
  prob_robust = 0.95,
  keep_el = TRUE,
  slide=1000, ## slide~400 is needed
  cores = 10,ld_method = "corr"
)
# LD_decay_emp <- readRDS("~/gitlab/formica_hybrid_old/data/ld_decay_corr.rds")
plot(ld_decay_sim)
plot(ld_decay_sim,type="chr",chr="Chr1")
plot(ld_decay_DI,type="chr",chr="Chr1")
plot(ld_decay_sim$decay_sum$a[-length(ld_decay_sim$decay_sum$a)],ld_decay_DI$decay_sum$a)
abline(0,1)
#ld_decay_sim$by_chr$Chr26$decay[,plot((start+end)/2,a,type="l")]
#ld_decay_DI$by_chr$Chr26$decay[,lines((start+end)/2,a,type="l",col="red")]

#dim(GTs_sim)
map[,table(DiagnosticIndex> -25)]
# plot(ld_decay)
# plot(ld_decay,type="chr",chr="ch26")
# ld_decay$decay_sum

ld_ws <- compute_ld_w(rho=0.95,ld_decay_sim)
names(ld_ws) <- map$marker
length(ld_ws)

plot(ld_ws,type="l")

p_0 <- c(rep(0.01,100),rep(0.5,40),rep(0.05,100))
image(cbind(sapply(1:10,function(x)rbinom(210,size = 1,p_0)),
      sapply(1:10,function(x)rbinom(210,size = 1,1-p_0))))


###### from old #######

rec_map <- fread("./data/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]


plot_ld_decay_tracks <- function(
    ld_decay,
    ld_w_095 = NULL,
    rec_dt = NULL,
    map_CeMLG = NULL,
    outliers = NULL,
    plot_ld_w = TRUE,
    plot_rec_rate = FALSE,
    highlight_outliers = TRUE,
    a_threshold = 0.0025,
    ncol = 5
) {
  
  library(data.table)
  library(ggplot2)
  
  ## Build LD-decay window table
  
  plot_dt <- rbindlist(lapply(names(ld_decay$by_chr), function(ch) {
    
    chr_obj <- ld_decay$by_chr[[ch]]
    decay <- copy(chr_obj$decay)
    
    decay[, Chr := ch]
    decay[, mid := rowMeans(.SD, na.rm = TRUE),
          .SDcols = c("start", "end")]
    decay[, chr_mean := chr_obj$decay_sum$a[1]]
    
    decay
  }), fill = TRUE)
  
  plot_dt[, Chr := factor(Chr, levels = names(ld_decay$by_chr))]
  
  ## Add median ld_w per window
  if (plot_ld_w && !is.null(ld_w_095)) {
    
    ld_dt <- data.table(
      snp = names(ld_w_095),
      ld_w = as.numeric(ld_w_095)
    )
    
    ld_dt[, c("Chr", "Pos") := tstrsplit(snp, ":", fixed = TRUE)]
    ld_dt[, Pos := as.numeric(Pos)]
    
    plot_dt[, ld_w_med := NA_real_]
    
    for (ch in names(ld_decay$by_chr)) {
      
      ld_ch <- ld_dt[Chr == ch]
      
      plot_dt[Chr == ch, ld_w_med := sapply(seq_len(.N), function(i) {
        median(
          ld_ch[Pos >= start[i] & Pos <= end[i], ld_w],
          na.rm = TRUE
        )
      })]
    }
    
    ## scale ld_w to comparable axis as a
    plot_dt[, scale_fac_ld := max(a, na.rm = TRUE) /
              max(ld_w_med, na.rm = TRUE),
            by = Chr]
    
    plot_dt[, ld_w_scaled := ld_w_med * scale_fac_ld]
  }
  
  ## Add recombination rate per LD-decay window
  if (plot_rec_rate && !is.null(rec_dt)) {
    
    rec_dt <- copy(rec_dt)
    
    if (!"Chr" %in% names(rec_dt)) {
      rec_dt[, Chr := paste0("Chr", sub("chromosome_", "", chr))]
    }
    
    setDT(rec_dt)
    
    plot_dt[, win_id := .I]
    
    rec_win <- rec_dt[
      plot_dt,
      on = .(Chr, pos >= start, pos <= end),
      allow.cartesian = TRUE
    ][
      ,
      .(
        rec_rate = mean(`cM/Mb`, na.rm = TRUE),
        rec_cM   = mean(cM, na.rm = TRUE)
      ),
      by = win_id
    ]
    
    plot_dt <- rec_win[plot_dt, on = "win_id"]
    
    plot_dt[, scale_fac_rec := max(a, na.rm = TRUE) /
              max(rec_rate, na.rm = TRUE),
            by = Chr]
    
    plot_dt[, rec_scaled := rec_rate * scale_fac_rec]
  }
  
  ## Base plot
  p <- ggplot(plot_dt, aes(mid / 1e6)) +
    geom_line(aes(y = a), linewidth = 0.6, col = "salmon") +
    facet_wrap(~ Chr, scales = "free", ncol = ncol) +
    labs(
      x = "Chromosome position (Mbp)",
      y = "Scaled tracks"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ## Add threshold
  if (!is.null(a_threshold)) {
    p <- p +
      geom_hline(
        yintercept = a_threshold,
        linetype = 2,
        linewidth = 0.5
      )
  }
  
  ## Add LD-cluster SNPs
  if (!is.null(map_CeMLG)) {
    
    map_CeMLG <- copy(map_CeMLG)
    
    p <- p +
      geom_point(
        data = map_CeMLG[color != "grey"],
        aes(x = Pos / 1e6, y = ld_w_095 / 100),
        inherit.aes = FALSE,
        size = 0.25,
        alpha = 0.5,
        col = "grey40"
      )
  }
  
  ## Add window median ld_w
  if (plot_ld_w && "ld_w_scaled" %in% names(plot_dt)) {
    p <- p +
      geom_line(
        aes(y = ld_w_scaled),
        linewidth = 0.4,
        col = "grey30",
        alpha = 0.8
      )
  }
  
  ## Add recombination rate
  if (plot_rec_rate && "rec_scaled" %in% names(plot_dt)) {
    p <- p +
      geom_line(
        aes(y = rec_scaled),
        linewidth = 0.5,
        col = "steelblue",
        alpha = 0.8
      )
  }
  
  ## Highlight outliers
  if (highlight_outliers && !is.null(outliers)) {
    
    outliers <- copy(outliers)
    
    if (!"mid" %in% names(outliers)) {
      outliers[, mid := rowMeans(.SD, na.rm = TRUE),
               .SDcols = c("start", "end")]
    }
    
    p <- p +
      geom_point(
        data = outliers,
        aes(x = mid / 1e6, y = a),
        inherit.aes = FALSE,
        size = 1,
        col = "red"
      )
  }
  
  p
}

map[,color:="grey"]
p <- plot_ld_decay_tracks(
  ld_decay = ld_decay,
  ld_w_095 = ld_ws[,"0.95"],
  rec_dt = rec_map,
  map_CeMLG = NULL,
  plot_ld_w = TRUE,
  plot_rec_rate = TRUE,
  outliers = NULL,
  highlight_outliers = FALSE
)
p

## ----------

ld_dt <- data.table(
  snp = names(ld_w_095),
  ld_w = as.numeric(ld_w_095)
)

ld_dt[, c("Chr", "Pos") := tstrsplit(snp, ":", fixed = TRUE)]
ld_dt[, Pos := as.numeric(Pos)]


rec_map <- copy(rec_map)

rec_map[, `:=`(
  d_bp = shift(pos, type = "lead") - pos,
  d_cM = shift(cM, type = "lead") - cM
), by = Chr]

rec_map[, rec_rate := d_cM / (d_bp / 1e6)]


rec_map[!is.na(rec_rate),mid := (pos + shift(pos, type = "lead")) / 2,by = Chr]

rec_int <- rec_map[!is.na(d_bp),.(Chr, pos, mid, rec_rate,`cM/Mb`)]

setkey(rec_int, Chr, pos)
setkey(ld_dt, Chr, Pos)

ld_rec <- rec_int[ld_dt, roll = TRUE]

rec_summary <- ld_rec[
  ,
  .(
    median_ldw = median(ld_w),
    n_snps = .N
  ),
  by = .(Chr, mid, rec_rate,`cM/Mb`)
]

par(mfcol=c(2,1))
rec_summary[,plot(median_ldw,rec_rate)]
rec_summary[,plot(median_ldw,`cM/Mb`)]

dt <- melt(rec_summary,measure.vars = c("cM/Mb","rec_rate"),variable.name = "Rec_type")


p1 <- ggplot(dt[],
             aes(value, median_ldw,col=log10(n_snps))) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(option="turbo")+
  geom_smooth(method = "loess",se=FALSE) +
  facet_grid(Rec_type~.) +
  theme_bw(base_size = 12)
p

p2 <-ggplot(dt[],
            aes(log10(n_snps),value,col=median_ldw)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c(option="turbo")+
  geom_smooth(method = "lm",se=FALSE) +
  facet_grid(Rec_type~.) +
  theme_bw(base_size = 12)

p1 | p2


rec_map[,plot(log10(d_bp), `cM/Mb`)]
ggplot(rec_map, aes(log10(d_bp), `cM/Mb`)) +
  geom_point()+
  theme_bw(base_size = 12)+
  geom_smooth(method = "lm",se=FALSE)

ggplot(dt) +
  geom_line(data=dt[Rec_type=="rec_rate"], aes(mid, value/2),col="grey")+
  geom_line(data=dt[Rec_type!="rec_rate"], aes(mid, value),col="black")+
  facet_wrap(Chr~.,scales="free_x") +
  theme_bw(base_size = 12) +
  ylim(0,200)

p1 <- ggplot(rec_summary, aes(mid,`cM/Mb`*median_ldw)) +
  geom_point()+
  facet_wrap(Chr~.,scales="free_x") +
  theme_bw(base_size = 12)

p2 <- ggplot(rec_summary, aes(mid,rec_rate*median_ldw)) +
  geom_point()+
  facet_wrap(Chr~.,scales="free_x") +
  theme_bw(base_size = 12)

p1 / p2

par(mfcol=c(2,1))
rec_summary[,hist(log10(1+rec_rate)*median_ldw,breaks=100,xlim=c(0,3))]
rec_summary[,hist(log10(1+`cM/Mb`),breaks=40,xlim=c(0,3))]

#rec_summary
rec_summary[,rec_q:=ecdf(rec_rate)(rec_rate)]
rec_summary[,ld_q:=ecdf(median_ldw)(median_ldw)]
rec_summary[,s:=rec_q * ld_q]
rec_summary[,quantile(s,0.99,na.rm=TRUE)]
rec_summary[,hist(s)]

ggplot(rec_summary,aes(mid, s))+
  geom_point()+
  facet_wrap(Chr~.,scales="free_x") +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0.75,col="salmon")

ld_rec[,rec_q:=ecdf(`cM/Mb`)(`cM/Mb`)]
#ld_rec[,rec_q:=ecdf(rec_rate)(rec_rate)]
ld_rec[,ld_q:=ecdf(ld_w)(ld_w)]
ld_rec[,s:=rec_q * ld_q]
ld_rec[,hist(s)]

ld_rec[,rec_rate_norm:=rec_rate/max(rec_rate,na.rm=TRUE)]

ggplot(ld_rec,aes(pos, ld_w*100))+
  geom_point(size=0.25,alpha=0.5,col="grey")+
  geom_line(data=dt[Rec_type=="rec_rate"], aes(mid, value),col="black")+
  facet_wrap(Chr~.,scales="free_x") +
  theme_bw(base_size = 12) #+
#geom_hline(yintercept = 0.75,col="salmon")

rec_sm <- copy(rec_map)

fit <- rec_sm[!is.na(rec_rate),smooth.spline(pos, cM, spar=0.5)]

rec_sm[!is.na(rec_rate),rec_rate_pred := predict(smooth.spline(pos, cM, spar=0.1),x = pos,deriv = 1)$y*1e6,by=Chr]

ggplot(ld_rec,aes(pos, ld_w*100))+
  geom_point(size=0.25,alpha=0.5,col="grey")+
  geom_line(data=rec_sm, aes(mid, `cM/Mb`),col="black",inherit.aes = FALSE)+
  facet_wrap(Chr~.) +
  theme_bw(base_size = 12) #+

ggplot(rec_sm,aes(rec_rate_pred, ld_w*100))+
  geom_point(size=0.25,alpha=0.5,col="grey")+
  facet_wrap(Chr~.,scales="free_x") +
  theme_bw(base_size = 12) #+


dt <- rec_summary[is.finite(`cM/Mb`) & !is.na(`cM/Mb`)]

dt[, log_rec := log10(1 + `cM/Mb`)]

dt[,plot(log_rec,median_ldw)]
#library(mgcv)
m1 <- gam(
  log_rec ~ s(median_ldw, k = 10),
  data = dt[Chr=="Chr1"],
  method = "REML"
)

summary(m1)

dt[Chr=="Chr1", pred_log_rec := predict(m1)]

plot(dt[Chr=="Chr1"]$log_rec, dt[Chr=="Chr1"]$pred_log_rec, use = "complete.obs")^2

rec_rate <- pred$y * 1e6
plot(rec_rate)

ggsave("tmp.png",height = 9,width = 8,units = "in",dpi = 300,plot = p)

setkey(rec_map, Chr, mid)

plot_dt_rec <- rec_map[
  plot_dt,
  on = .(Chr, pos=mid),
  roll = "nearest"
]

dt <- melt(plot_dt_rec,measure.vars = c("cM/Mb","rec_rate"),variable.name = "Rec_type")
ggplot(dt,
       aes(value, ld_w_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  facet_grid(Rec_type~.)


plot_dt[, win_id := .I]

rec_window <- rec_map[
  plot_dt,
  on = .(Chr,
         pos >= start,
         pos <= end),
  allow.cartesian = TRUE
]

rec_window <- rec_window[
  ,
  .(
    rec_rate = mean(`cM/Mb`, na.rm = TRUE),
    rec_cM   = mean(cM, na.rm = TRUE)
  ),
  by = win_id
]

plot_dt_rec <- rec_window[plot_dt, on = "win_id"]

ggplot(plot_dt_rec,
       aes(rec_cM, ld_w_med)) +
  geom_point(alpha = 0.5)

library(quantreg)
library(splines)
ok <- complete.cases(plot_dt_rec[, .(ld_w_med, rec_rate)])

fit95 <- rq(
  ld_w_med ~ splines::bs(rec_rate, df = 7),
  tau = 0.95,
  data = plot_dt_rec[ok]
)

plot_dt_rec[, ld95 := NA_real_]
plot_dt_rec[ok, ld95 := predict(fit95, newdata = plot_dt_rec[ok])]

plot_dt_rec[, excess_ld := ld_w_med - ld95]

candidates <- plot_dt_rec[
  excess_ld > quantile(excess_ld, 0.99, na.rm = TRUE)
]



ggplot(plot_dt_rec, aes(rec_rate, ld_w_med)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(
    data = candidates,
    color = "red",
    size = 0.8
  ) +
  theme_bw(base_size = 12) +
  geom_line(
    data = plot_dt_rec[ok][order(rec_rate)],
    aes(y = ld95),
    linewidth = 1
  )

##
ok <- complete.cases(plot_dt_rec[, .(a, rec_rate)])

fit95 <- rq(
  a ~ splines::bs(rec_rate, df = 7),
  tau = 0.95,
  data = plot_dt_rec[ok]
)

plot_dt_rec[, ld95 := NA_real_]
plot_dt_rec[ok, ld95 := predict(fit95, newdata = plot_dt_rec[ok])]

plot_dt_rec[, excess_ld := ld_w_med - ld95]

candidates <- plot_dt_rec[
  excess_ld > quantile(excess_ld, 0.99, na.rm = TRUE)
]


p1 <- ggplot(plot_dt_rec, aes(a, ld_w_med)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_bw(base_size = 12)

p2 <- ggplot(plot_dt_rec, aes(rec_rate,a)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_bw(base_size = 12)

p3<- ggplot(plot_dt_rec, aes(rec_rate,ld_w_med)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_bw(base_size = 12) +
  geom_smooth(se=FALSE)

p4<- ggplot(plot_dt_rec, aes(rec_rate,sample(ld_w_med))) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_bw(base_size = 12)+
  geom_smooth(se=FALSE)


p1 | p2 | p3

p3+ggtitle("Empirical data") | p4+ggtitle("ld_w permuted")


###### even older #######
ld_w_0.95 <- ld_ws[,"0.95"]
names(ld_w_0.95) <- map_1mb$marker

plot_dt <- rbindlist(
  lapply(names(ld_decay$by_chr), function(ch) {
    
    chr_obj <- ld_decay$by_chr[[ch]]
    
    decay <- copy(chr_obj$decay)
    
    decay[, mid := rowMeans(.SD, na.rm = TRUE),
          .SDcols = c("start", "end")]
    
    decay[, chr := ch]
    
    decay[, genome_mean := chr_obj$decay_sum$a[1]]
    
    decay
  })
)

## optional: chromosome ordering
plot_dt[, chr := factor(
  chr,
  levels = paste0("Chr", 1:25)
)]

## plot
# ggplot(plot_dt, aes(mid, a)) +
# 
#   geom_line(linewidth = 0.4) +
# 
#   geom_hline(
#     aes(yintercept = genome_mean),
#     linetype = 2,
#     linewidth = 0.5
#   ) +
# 
#   facet_wrap(~ chr, scales = "free_x", ncol = 5) +
# 
#   labs(
#     x = "Chromosome position",
#     y = "Decay rate (a)"
#   ) +
# 
#   theme_bw(base_size = 10) +
# 
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_text(face = "bold"),
#     panel.grid.minor = element_blank(),
#     panel.spacing = unit(0.8, "lines")
#   )
#



## with ld_w

ld_dt <- data.table(
  snp = names(ld_w_0.95),
  ld_w = as.numeric(ld_w_0.95)
)

ld_dt[, c("Chr", "Pos") := tstrsplit(snp, ":", fixed = TRUE)]
ld_dt[, Pos := as.numeric(Pos)]

## add local ld_w median to each decay window
plot_dt <- rbindlist(lapply(names(ld_decay$by_chr), function(ch) {
  
  chr_obj <- ld_decay$by_chr[[ch]]
  decay <- copy(chr_obj$decay)
  
  decay[, Chr := ch]
  decay[, mid := rowMeans(.SD, na.rm = TRUE),
        .SDcols = c("start", "end")]
  decay[, chr_mean := chr_obj$decay_sum$a[1]]
  
  ld_ch <- ld_dt[Chr == ch]
  
  decay[, ld_w_med := sapply(seq_len(.N), function(i) {
    median(ld_ch[Pos >= start[i] & Pos <= end[i], ld_w],
           na.rm = TRUE)
  })]
  
  decay
}), fill = TRUE)

plot_dt[, Chr := factor(Chr, levels = names(ld_decay$by_chr))]

plot_dt[, scale_fac_chr := max(a, na.rm = TRUE) / max(ld_w_med, na.rm = TRUE),by = Chr]

plot_dt[, ld_w_scaled := ld_w_med * scale_fac_chr]

#plot_dt[is.na(Chr)]



p_LD_decay_local_LD <- ggplot(plot_dt, aes(mid/1e6)) +
  #geom_point(data=map_CeMLG[color=="grey"],aes(Pos/1e6,ld_w_0.95/100),size=0.1,alpha=0.5)+
  #geom_point(data=map_CeMLG[color!="grey"],aes(Pos/1e6,ld_w_0.95/100),size=0.25,alpha=0.5)+
  geom_line(aes(y = a), linewidth = 0.6,col="steelblue") +
  geom_hline(aes(yintercept = 0.0025),
             linetype = 2, linewidth = 0.6) +
  facet_wrap(~ Chr, scales = "free", ncol = 5) +
  labs(
    x = "Chromosome position (Mbp)",
    y = expression("Decay rate (a) | "*ld["w,"*rho*"=0.95"]/100)
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
p_LD_decay_local_LD

ggplot(plot_dt, aes(a,ld_w_med)) +
  geom_point() +
  theme_bw(base_size = 22) +
  geom_smooth(se=F) +
  ylab(expression(ld["w,"*rho*"=0.95"]~"(in windows)"))+
  theme(panel.grid.minor = element_blank()) +
  ggtitle(expression("Local LD "*(ld["w"])*" vs. decay rate"))



