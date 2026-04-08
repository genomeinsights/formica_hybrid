library(LDna)
library(Matrix)
library(wesanderson)
library(ggplot2)
library(data.table)
library(SNPRelate)
library(LDscnR)
library(circlize)

col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
                "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
                "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
                "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")

data <- readRDS("./data/filtered_DIEM_1mb.rds")
map <- data$map
GTs <- data$GTs
sample_info <- data$sample_info
rm(data)
gc()

ld_decay <- readRDS("./data/ld_decay_corr.rds")

gds <- create_gds_from_geno(geno=GTs, map, "gds_formicia")

map <- readRDS("LDscnR_data.rds")$map

map[,table(DiagnosticIndex<(-15),C_joint>0.2)]
map[DiagnosticIndex>(-15),table(C_joint>0.2)]

#ld_w_99 <- compute_ld_w(ld_decay,rho = 0.99)

#plot(ld_w_99)
#ch ="Chr2"
#27/16
# rho_d = 0.95
# rho_ld = 0.9
# par(mfcol=c(4,4))
# par(mar=c(1,1,1,1))
#
# MAP = map[ld_w_95>0.1 | C_joint>0.2 ]

find_hapl_blocks <- function(ld_decay,map,SNPs,rho_ld=0.95,rho_d=0.95,ld_th=NULL,d_th=NULL,col_vector=NULL){
  #ch = "Chr10"
  if(is.null(col_vector)){
    col_vector <- c("#B2DF8A", "#FFD92F", "firebrick", "#33A02C", "#7FC97F", "#CAB2D6",
                    "#FB8072", "grey30", "#E6AB02", "#FDC086", "steelblue", "#1F78B4",
                    "#FB9A99", "#1B9E77", "#BC80BD", "#E31A1C", "#7570B3", "#A6761D",
                    "#A6CEE3", "salmon", "#FFFF33", "forestgreen", "#FDCDAC", "#BF5B17",
                    "#A6761D", "#FBB4AE", "#4DAF4A", "#B3E2CD", "#FDDAEC", "#BEBADA",
                    "#FFF2AE", "#1F78B4", "#66C2A5", "#F0027F", "#E6AB02", "#E78AC3",
                    "#FF7F00", "#8DA0CB", "#6A3D9A", "#B15928", "#E41A1C")
  }

  map_with_hb <- rbindlist(lapply(names(ld_decay$by_chr), function(ch){
    message(ch)
    mp <- copy(map[Chr==ch])

    chr_obj <- ld_decay$by_chr[[ch]]
    chr_obj$el <- fread(chr_obj$el,showProgress = FALSE)

    if(is.null(ld_th) & is.null(d_th)){
      a_chr <- ld_decay$decay_sum[Chr==ch,a_pred]
      b_chr <- ld_decay$decay_sum[Chr==ch,b]
      c_chr <- ld_decay$decay_sum[Chr==ch,c_pred]

      d_th  <- d_from_rho(a_chr, rho = rho_d)
      ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_ld)
      ed <- chr_obj$el[r2>ld_th & d<d_th & (SNP1 %in% SNPs | SNP2 %in% SNPs),.(SNP1,SNP2)]
    }else{
      ed <- chr_obj$el[r2>ld_th & d<d_th & (SNP1 %in% SNPs | SNP2 %in% SNPs),.(SNP1,SNP2)]
    }



    #print(dim(ed))
    g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
    comps <- igraph::components(g)

    ors_chr <- split(names(comps$membership), comps$membership)
    ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= 10]

    hap_bl <- data.table(HB_id=rep(paste0(ch,"_",1:length(ors_chr)),sapply(ors_chr,length)),marker=unlist(ors_chr),col=rep(rep(sample(col_vector),100)[1:length(ors_chr)],sapply(ors_chr,length)))
    mp[,HB_id:=hap_bl$HB_id[match(marker,hap_bl$marker)]]
    mp[,HB_col:=hap_bl$col[match(marker,hap_bl$marker)]]
    mp[is.na(HB_id),HB_col:="grey"]

    return(mp)
  }))
  return(map_with_hb)
}

map_HB_095_090 <- find_hapl_blocks(ld_decay,map=map,SNPs=map[ld_w_99>0.1 | C_joint>0.2 ,marker],rho_ld=0.9,rho_d=0.95)
#mp[!is.na(OR_id),plot(Pos,ld_w_95,col=OR_col,cex=0.25,xlim=range(mp[,Pos]))]
map_HB_095_090[,Chr := factor(Chr, levels = unique(Chr)[order(as.numeric(gsub("Chr","",unique(Chr))))])]

map[ ,plot(ld_w_99,cex=0.2)]

map_HB_095_090[]
ggplot(map_HB_095_090[!is.na(HB_id)], aes(Pos/1e6,ld_w_99,col=HB_col)) +
  geom_point(size=0.25) +
  scale_color_identity() +
  facet_wrap(Chr~.,scale="free_x",nrow=4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank()) +
  geom_hline(yintercept = 0.2,linetype=2)

map_HB_099_095[,]
layout <- prep_manhattan(map_HB_099_095[C_joint>0.2 ])

p_PC1 <- plot_manhattan_gg(layout, y_vars = c("C_joint","C_PC1","C_PC2"), y_labels = c("C_joint","C_PC1","C_PC2"),
                           thresholds = c(0.2,0.2,0.2), col_var = "HB_col",
                           point_size = 1, ncol = 1)

ggplot(map_HB_099_095[!is.na(HB_id)], aes(Pos/1e6,ld_w_99,col=HB_col)) +
  geom_point(size=0.25) +
  scale_color_identity() +
  facet_wrap(Chr~.,scale="free_x",nrow=4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank()) +
  geom_hline(yintercept = 0.2,linetype=2)

# map_ORs[,indx := .I]
# keep <- unlist(map_ORs[,.(indx[which.max(Pos)],indx[which.max(Pos)]),by=Chr][,.(V1,V1)])
# dt <- rbind(map_ORs[!is.na(OR_id)],map_ORs[keep])
#
# ggplot(dt,aes(Pos/1e6, ld_w_99,col=OR_col))+
#   geom_point(size=0.25) +
#   facet_wrap(Chr~.,scale="free_x",nrow=4)+
#   scale_color_identity() +
#   theme_minimal()

saveRDS(map_HBs,"map_HBs.rds")
map_HBs <- map_HB_095_09
core_SNPs <- map_HBs[!is.na(HB_id),.(marker=marker[which.max(ld_w_95*pmax(F_BF_PC1,   F_BF_PC2,na.rm = TRUE))],ld_w=max(ld_w_95),N=.N),by=HB_id]

map_HBs[,core := marker %in% core_SNPs[,marker]]

table(map_HBs[C_PC1>0.2,unique(HB_id)] %in% map_HBs[C_PC2>0.2,unique(HB_id)])
table( map_HBs[C_PC2>0.2,unique(HB_id)] %in% map_HBs[C_PC1>0.2,unique(HB_id)])

map_HBs[C_PC1>0.2,unique(HB_id)][map_HBs[C_PC1>0.2,unique(HB_id)] %in% map_HBs[C_PC2>0.2,unique(HB_id)]]

map_HBs[C_PC1>0.2,table(HB_id)]
map_HBs[C_PC2>0.2,table(HB_id)]

do_these <- split(map_HBs$marker,map_HBs$HB_id)
hb <- do_these[[1]]
PCAs <- mclapply(do_these,function(hb){
  pca <- snpgdsPCA(gds,snp.id = hb,autosome.only = FALSE,verbose = FALSE)
  pca$eigenvect[,1]
},mc.cores=8)

PCAs <- do.call(cbind,PCAs)

heatmap(PCAs,scale="none")

cmat <- cor(PCAs)^2
cmat[upper.tri(cmat,diag = TRUE)] <- NA
#image(cmat)
dimnames(cmat) <- list(colnames(PCAs),colnames(PCAs))
#cmat[cmat<0.5] <- 0
ldna <- LDnaRaw(cmat)
cls_10 <- extractBranches(ldna,min.edges = 1)
#dev.off()
idx <- which(lower.tri(cmat), arr.ind = TRUE)

el <- data.table(
  cl1 = rownames(cmat)[idx[, 1]],
  cl2 = colnames(cmat)[idx[, 2]],
  r2   = cmat[idx],
  stringsAsFactors = FALSE
)
el$SNP1=map_HBs[which(core)][match(el$cl1,HB_id),marker]
el$SNP2=map_HBs[which(core)][match(el$cl2,HB_id,nomatch = 0),marker]
el$Chr1=map_HBs[which(core)][match(el$cl1,HB_id,nomatch = 0),Chr]
el$Chr2=map_HBs[which(core)][match(el$cl2,HB_id,nomatch = 0),Chr]
el$pos1=map_HBs[which(core)][match(el$cl1,HB_id,nomatch = 0),Pos]
el$pos2=map_HBs[which(core)][match(el$cl2,HB_id,nomatch = 0),Pos]
el[Chr1==Chr2,d:=abs(pos1-pos2)]
el[Chr1!=Chr2,d:=1e9]

#rbind(el[Chr1==Chr2 & d>1e6],el[Chr1!=Chr2 ])[r2>0.8]

# optional: remove missing or zero values
#el <- el[!is.na(el$r2), ]

tmp <- rbind(el[Chr1==Chr2 & d>1e6],el[Chr1!=Chr2 ])[r2>0.3][,unique(c(cl1,cl2))]
outl_HBs <- map_HBs[C_joint>0.2 & HB_id %in% tmp,unique(HB_id) ]



map_HBs[,OR_col:=HB_col]
map_HBs[,OR_id:=HB_id]
png("./Figures/circ_manh_out_095_09_03_5e6.png",height = 6,width = 6,res = 300,units = "in")
plot_ld_outlier_circos(
  el = el[(cl1 %in% outl_HBs | cl2 %in% outl_HBs) & r2>0.3 & d>5e6],
  #el = rbind(el[Chr1==Chr2 & d>5e6],el[Chr1!=Chr2 ])[r2>0.5],
  map = map_HBs,
  map_manh = map_HBs,
  r2_th = 0,
  only_inter = FALSE,
  add_manhattan = TRUE,
  manhattan_stat_col = c("ld_w_95","C_joint"),
  manhattan_transform = identity,
  manhattan_type = "points",
  pad_bp = 5e4,
  manhattan_color_col = "HB_col",
  manhattan_show_labels = FALSE,
  highlight_linked_or = FALSE,manhattan_cex = 0.15,manhattan_track_height = 0.05,
  show_legend = TRUE, outlier_col = "steelblue",
  linked_outlier_col = "firebrick4", link_alpha_range = c(1, 1),  link_lwd_range = c(2, 2),
  legend_x_offset = -0.01,
  legend_y_bottom = 0.35,
  legend_y_top = 0.75
)
dev.off()

connected_clusters <- el[(cl1 %in% outl_HBs | cl2 %in% outl_HBs) & r2>0.3 & d>5e6][,unique(c(cl1,cl2))]
keep <- which(colnames(PCAs) %in% connected_clusters)
PCA_connected <- PCAs[,keep]
rownames(PCA_connected) <- sample_info$Sample_ID
plot(cor(PCA_connected)[1,])

#plot(rowMeans(cor(PCA_connected)))
PCA_connected[,rowMeans(cor(PCA_connected))<0] <- -1*PCA_connected[,rowMeans(cor(PCA_connected))<0]
heatmap(PCA_connected,scale="none",cexRow=0.5)

colnames(PCA_connected)
Chr_un <- map_HBs[HB_id %in% colnames(PCA_connected),unique(Chr)]

cols <- col_vector[1:14][match(map_HBs[match(colnames(PCA_connected),map_HBs[,HB_id]),Chr],Chr_un)]
cols_pop <- col_vector[1:20][match(sample_info$Population,unique(sample_info$Population))]
png("./Figures/heatmap_PC1_outl_095_09_03_5e6.png",width = 8,height = 8,units = "in",res = 600)
heatmap(PCA_connected,scale="none",ColSideColors=cols,RowSideColors=cols_pop)
dev.off()
ncol(gts_connected_outliers)


keep_SNPs <- map_HBs[,which(HB_id %in% connected_clusters)]
gts_connected_outliers <- GTs[,keep_SNPs]
#image(cor(gts_connected_outliers[,1:1000],use="pair"))
polarity <- rowMeans(cor(gts_connected_outliers,use="pair"))
plot(polarity)
gts_connected_outliers[,polarity<0] <- 2-gts_connected_outliers[,polarity<0]

cols <- col_vector[1:27][match(map_HBs[marker %in% colnames(gts_connected_outliers),Chr],unique(map_HBs[marker %in% colnames(gts_connected_outliers),Chr]))]
cols_pop <- col_vector[1:20][match(sample_info$Population,unique(sample_info$Population))]
png("./Figures/heatmap_SNP_outl_095_09_03_5e6.png",width = 8,height = 8,units = "in",res = 600)
heatmap(gts_connected_outliers,scale="none",ColSideColors=cols,RowSideColors=cols_pop)
dev.off()
ncol(gts_connected_outliers)


length(map_HBs[marker %in% colnames(gts_connected_outliers),Chr])


png("./Figures/heatmap_SNP_outl_095_09_03_5e6_RowvNA.png",width = 8,height = 8,units = "in",res = 600)
heatmap(gts_connected_outliers,scale="none",cexRow=0.5,Rowv=NA)
dev.off()

png("./Figures/heatmap_SNP_outl_095_09_03_5e6_RowvNA_ColvNA.png",width = 8,height = 8,units = "in",res = 600)
heatmap(gts_connected_outliers,scale="none",cexRow=0.5,Colv=NA)
dev.off()

image(t(gts_connected_outliers))


####
map_HBs <- map_HB_095_09
core_SNPs <- map_HBs[!is.na(HB_id),.(marker=marker[which.max(ld_w_95*F_BF_PC2)],ld_w=max(ld_w_95),N=.N),by=HB_id]

map_HBs[,core := marker %in% core_SNPs[,marker]]


el <- get_el(gds, slide_win_ld=-1,idx = map_HBs[,which(core)],by_chr = FALSE,method = "corr")

table(map_HBs[C_PC1>0.2,unique(HB_id)] %in% map_HBs[C_PC2>0.2,unique(HB_id)])

connected_SNPs <- rbind(el[Chr1==Chr2 & d>5e6],el[Chr1!=Chr2 ])[r2>0.5][,unique(c(SNP1,SNP2))]


snps <- unique(c(el$SNP1, el$SNP2))
snps <- sort(snps)
n <- length(snps)

i <- match(el$SNP1, snps)
j <- match(el$SNP2, snps)

mat <- sparseMatrix(i = c(i, j),
                           j = c(j, i),
                           x = c(el$r2, el$r2),
                           dims = c(n, n),
                           dimnames = list(snps, snps))

mat <- as.matrix(tril(mat))
mat[upper.tri(mat,diag = TRUE)] <- NA




keep <- !is.na(colnames(mat))
mat <- mat[keep,keep]
#mat[1:10,1:10]
#heatmap(mat)
#dev.off()
ldna <- LDnaRaw(mat)
cls_10 <- extractBranches(ldna,min.edges = 5)
length(cls_10)

png("./Figures/circ_manh.png",height = 6,width = 6,res = 300,units = "in")
plot_ld_outlier_circos(
  el = rbind(el[Chr1==Chr2 & d>5e6],el[Chr1!=Chr2 ])[r2>0.3],
  map = map_HBs,
  map_manh = map_HBs,
  r2_th = 0,
  only_inter = FALSE,
  add_manhattan = TRUE,
  manhattan_stat_col = c("ld_w_95","C_joint"),
  manhattan_transform = identity,
  manhattan_type = "points",
  pad_bp = 5e4,
  manhattan_color_col = "HB_col",
  manhattan_show_labels = FALSE,
  highlight_linked_or = FALSE,manhattan_cex = 0.15,manhattan_track_height = 0.05,
  show_legend = TRUE, outlier_col = "steelblue",
  linked_outlier_col = "firebrick4", link_alpha_range = c(1, 1),  link_lwd_range = c(2, 2),
  legend_x_offset = -0.01,
  legend_y_bottom = 0.35,
  legend_y_top = 0.75
)
dev.off()


ed <-     rbind(el[Chr1==Chr2 & d>5e6],el[Chr1!=Chr2 ])[r2>0.5,.(SNP1,SNP2,r2)]
g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
comps <- igraph::components(g)
ors_chr <- split(names(comps$membership), comps$membership)

unlist(ors_chr)

el <- get_el(gds, slide_win_ld=-1,idx = map_ORs[,which(marker %in% unlist(ors_chr))],by_chr = FALSE,method = "corr")
el[,hist(r2)]

heatmap(GTs[,map_ORs[,which(marker %in% unlist(ors_chr))]],scale="none")


tmp <- cor(GTs[,map_ORs[,which(marker %in% unlist(ors_chr[2]))]],use="pair")^2
diag(tmp) <- NA
hist(tmp)

nodes <- unique(rbind(
  el[, .(idx = Var1, snp = SNP1)],
  el[, .(idx = Var2, snp = SNP2)]
))

setorder(nodes, idx)

n <- max(nodes$idx)
mat <- matrix(0, nrow = n, ncol = n,
              dimnames = list(nodes$snp[match(seq_len(n), nodes$idx)],
                              nodes$snp[match(seq_len(n), nodes$idx)]))

#ld_decay_noora
mat[cbind(el$Var1, el$Var2)] <- el$r2
diag(mat) <- NA
mat[upper.tri(mat)] <- NA
dimnames(mat) <- list(gsub(":","_",colnames(mat)),gsub(":","_",colnames(mat)))
keep <- !is.na(colnames(mat))
mat <- mat[keep,keep]
#image(mat)

ldna <- LDnaRaw(mat)
cls_10 <- extractBranches(ldna,min.edges = 1)
dim(GTs)
GTs
map_ORs[,LG:=as.numeric(gsub("Chr","",map_ORs$Chr))]
setorder(map_ORs,LG,Pos)
heatmap(GTs[,map_ORs[,which(marker %in% unlist(ors_chr))]],scale="none")
dev.off()

hist(cor(GTs[,map_ORs[,which(marker %in% unlist(ors_chr))]],use = "pair")^2)

image(t(GTs[,map_ORs[,which(marker %in% ors_chr[[2]])]]))
