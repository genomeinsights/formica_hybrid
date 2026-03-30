library(LDna)
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
get_el

data <- readRDS("./data/filtered_DIEM_1mb.rds")
map <- data$map
GTs <- data$GTs
sample_info <- data$sample_info
rm(data)
gc()

ld_decay <- readRDS("./data/ld_decay_corr.rds")

gds <- create_gds_from_geno(geno=GTs, map, "gds_formicia")

map[,ld_w_99:=compute_ld_w(ld_decay,rho = 0.99)]

#ch ="Chr2"
rho_d = 0.95
rho_ld = 0.75
par(mfcol=c(3,4))
par(mar=c(1,1,1,1))
map_ORs <- rbindlist(lapply(names(ld_decay$by_chr), function(ch){
  mp <- map[Chr==ch]
  a_chr <- ld_decay$decay_sum[Chr==ch,a_pred]
  b_chr <- ld_decay$decay_sum[Chr==ch,b]
  c_chr <- ld_decay$decay_sum[Chr==ch,c_pred]

  d_th  <- d_from_rho(a_chr, rho = rho_d)
  ld_th <- ld_from_rho(b_chr, c_chr, rho = rho_ld)

  #d_window <- d_from_rho(a, rho_d)

  chr_obj <- ld_decay$by_chr[[ch]]
  chr_obj$el <- fread(chr_obj$el,showProgress = FALSE)

  keep <- mp[ld_w_99>0.05,marker]
  ed <- chr_obj$el[r2>ld_th & d<d_th & SNP %in% keep,.(SNP1=paste(ch,pos1,sep=":"),SNP2=paste(ch,pos2,sep=":"))]
  g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
  comps <- igraph::components(g)

  ors_chr <- split(names(comps$membership), comps$membership)
  ors_chr <- ors_chr[vapply(ors_chr, length, integer(1)) >= 10]

  ORs <- data.table(OR_id=rep(paste0("OR_",ch,1:length(ors_chr)),sapply(ors_chr,length)),marker=unlist(ors_chr),col=rep(rep(col_vector,100)[1:length(ors_chr)],sapply(ors_chr,length)))
  mp[,OR_id:=ORs$OR_id[match(marker,ORs$marker)]]
  mp[,OR_col:=ORs$col[match(marker,ORs$marker)]]
  mp[is.na(OR_id),OR_col:="black"]

  mp[!is.na(OR_id),plot(Pos,ld_w_99,col=OR_col,cex=0.25,xlim=range(mp[,Pos]))]
  return(mp)

}))
# map_ORs[,indx := .I]
# keep <- unlist(map_ORs[,.(indx[which.max(Pos)],indx[which.max(Pos)]),by=Chr][,.(V1,V1)])
# dt <- rbind(map_ORs[!is.na(OR_id)],map_ORs[keep])
#
# ggplot(dt,aes(Pos/1e6, ld_w_99,col=OR_col))+
#   geom_point(size=0.25) +
#   facet_wrap(Chr~.,scale="free_x",nrow=4)+
#   scale_color_identity() +
#   theme_minimal()


core_SNPs <- map_ORs[!is.na(OR_id),.(marker=marker[which.max(ld_w_99)],ld_w=max(ld_w_99),N=.N),by="OR_id"]

map_ORs[,core := marker %in% core_SNPs[,marker]]

el <- get_el(gds, slide_win_ld=-1,idx = map_ORs[,which(core)],by_chr = FALSE,symmetric = FALSE,edge_symmetry = FALSE)

el[Chr1!=Chr2 & r2>0.8]
el[Chr1==Chr2 & r2>0.8 & d>1e6]

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

length(cls_10)

png("./circ_manh.png",height = 6,width = 6,res = 300,units = "in")
plot_ld_outlier_circos(
  el = rbind(el[Chr1==Chr2 & d>5e6],el[Chr1!=Chr2 ])[r2>0.8],
  map = map_ORs,
  map_manh = map_ORs,
  r2_th = 0,
  only_inter = FALSE,
  add_manhattan = TRUE,
  manhattan_stat_col = c("ld_w_99"),
  manhattan_transform = identity,
  manhattan_type = "points",
  pad_bp = 5e4,
  manhattan_color_col = "OR_col",
  manhattan_show_labels = FALSE,
  highlight_linked_or = FALSE,manhattan_cex = 0.15,manhattan_track_height = 0.05,
  show_legend = TRUE, outlier_col = "steelblue",
  linked_outlier_col = "firebrick4", link_alpha_range = c(1, 1),  link_lwd_range = c(2, 2),
  legend_x_offset = -0.01,
  legend_y_bottom = 0.35,
  legend_y_top = 0.75
)
dev.off()


ed <-     rbind(el[Chr1==Chr2 & d>5e6],el[Chr1!=Chr2 ])[r2>0.8,.(SNP1,SNP2,r2)]
g     <- igraph::graph_from_data_frame(ed, directed = FALSE)
comps <- igraph::components(g)
ors_chr <- split(names(comps$membership), comps$membership)

unlist(ors_chr)

el <- get_el(gds, slide_win_ld=-1,idx = map_ORs[,which(marker %in% unlist(ors_chr))],by_chr = FALSE,symmetric = FALSE,edge_symmetry = FALSE)
el[,hist(r2)]
heatmap(cor(GTs[,map_ORs[,which(marker %in% unlist(ors_chr[2]))]],use="pair")^2,scale="none")

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
