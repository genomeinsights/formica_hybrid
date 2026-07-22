library(data.table)
library(SNPRelate)
devtools::load_all("~/gitlab/LDscnR") 

if(!file.exists("./data/parsed_vcf.rds")){
  GT <- fread(
    "./data/large.012",
    header = FALSE,
    na.strings = "-1"
  )
  #dim(GT_ind)
  # Remove VCFtools' individual-number column
  GT[, V1 := NULL]
  
  samples <- fread(
    "./data/large.012.indv",
    header = FALSE
  )[[1]]
  
  map <- fread(
    "./data/large.012.pos",
    header = FALSE,
    col.names = c("chr", "pos")
  )
  table(map$chr)
  map <- data.table(
    Chr = gsub("chromosome_","Chr",map$chr),
    Pos = as.integer(map$pos)
  )
  
  dim(GT)
  
  map[,marker:=paste(Chr,Pos,sep=":")]
  colnames(GT) <- map$marker
  rownames(GT) <- samples
  saveRDS(list(GTs=GT,map=map),"./data/parsed_vcf.rds")  
  
}else{
  data <- readRDS("./data/parsed_vcf.rds")
  GTs <- data$GTs
  map <- data$map
  rm(data)
}

gds <- create_gds_from_geno(geno=GTs, map, "gds_vcf")

map[,maf:= snpgdsSNPRateFreq(gds)$MinorFreq]
# on.exit(rm(gds))
# on.exit(unlink(gds))
ld_decay_vcf <- compute_LD_decay(
  gds,
  el_data_folder = "./EL_vcf/",
  keep_el = TRUE,
  slide=400, 
  cores = 10,ld_method = "corr"
)

plot(ld_decay_vcf,type="chr",chr="Chr26")
ld_w_095 <- as.vector(compute_ld_w(0.95,ld_decay = ld_decay_vcf))
el <- fread(ld_decay_vcf$by_chr$Chr1$el)
chr_obj <- ld_decay_vcf$by_chr$Chr1
a <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,a_pred]
b <- ld_decay$decay_sum[Chr==chr_obj$decay_sum$Chr,b]

d_window <- d_from_rho(a, rho=0.95)

if(is.null(chr_obj$el)) stop("No edge list present")

if(is.character(chr_obj$el)) chr_obj$el <- fread(chr_obj$el,showProgress = FALSE)

#make symmetric
chr_obj$el <- data.table::rbindlist(list(
  chr_obj$el[, .(SNP = SNP1, pos = pos1, pos_other = pos2, r2, d)],
  chr_obj$el[, .(SNP = SNP2, pos = pos2, pos_other = pos1, r2, d)]
))

ld_w <- chr_obj$el[d<d_window,.(r2_median=median(r2)),by=SNP]

plot(ld_w[match(chr_obj$snp_ids,ld_w$SNP),r2_median])

length(ld_decay_vcf$by_chr$Chr1$el)
table(ld_decay_vcf$by_chr$Chr1$snp_ids==map[Chr=="Chr1",marker])


names(ld_w_095) <- map$marker
plot(ld_w_095)

plot_dt <- rbindlist(
  lapply(names(ld_decay$by_chr), function(ch) {
    
    chr_obj <- ld_decay$by_chr[[ch]]
    
    decay <- copy(chr_obj$decay)
    
    decay[, mid := rowMeans(.SD, na.rm = TRUE),.SDcols = c("start", "end")]
    
    decay[, chr := ch]
    
    decay[, genome_mean := chr_obj$decay_sum$a[1]]
    
    decay
  })
)

## optional: chromosome ordering
plot_dt[, chr := factor(
  chr,
  levels = paste0("Chr", 1:27)
)]
plot_dt <- plot_dt[is.na(a)]
plot_dt[is.na(a),a:=0]

#plot_dt[,cor.test(a,ld_w_med)]
## plot
#plot_dt[,a]
ggplot(plot_dt, aes(mid, a)) +

  geom_line(linewidth = 0.4) +

  geom_hline(
    aes(yintercept = genome_mean),
    linetype = 2,
    linewidth = 0.5
  ) +

  facet_wrap(~ chr, scales = "free_x", ncol = 5) +

  labs(
    x = "Chromosome position",
    y = "Decay rate (a)"
  ) +

  theme_bw(base_size = 10) +

  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )
#



## with ld_w

ld_dt <- data.table(
  snp = names(ld_w_095),
  ld_w = as.numeric(ld_w_095)
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
  geom_point(data=map,aes(Pos/1e6,ld_w_095/100),size=0.1,alpha=0.5)+
  #geom_point(data=map_CeMLG[color!="grey"],aes(Pos/1e6,ld_w_095/100),size=0.25,alpha=0.5,col="grey40")+
  geom_line(aes(y = a), linewidth = 0.6,col="salmon") +
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
  ylab(expression(ld["w,"*rho*"=0.95"]~"(in windows)"))+
  theme(panel.grid.minor = element_blank()) +
  ggtitle(expression("Local LD "*(ld["w"])*" vs. decay rate"))

## recombinatin map -------------------------------------------------------
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

p <- plot_ld_decay_tracks(
  ld_decay = ld_decay,
  ld_w_095 = ld_w_095,
  rec_dt = rec_map,
  map_CeMLG = map_CeMLG,
  plot_ld_w = TRUE,
  plot_rec_rate = TRUE,
  outliers = candidates,
  highlight_outliers = TRUE
)
p

## ----------

ld_dt <- data.table(
  snp = names(ld_w_095),
  ld_w = as.numeric(ld_w_095)
)

ld_dt[, c("Chr", "Pos") := tstrsplit(snp, ":", fixed = TRUE)]
ld_dt[, Pos := as.numeric(Pos)]

rec_map <- fread("./rec_map/Frufa_DTOL_PR.ref_genome.recmap")
rec_map[, Chr := paste0("Chr", sub("chromosome_", "", chr))]

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



