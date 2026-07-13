library(vcfR)
library(data.table)

VCF <- read.vcfR("./data/Formica_hybrids_20pops.vs.Frufa_DTOL_PR.ref_genome.unmasked.clipped.chromosomes.snpqual_readsupport.fixedHeader.indDP.missingness.biallSNPs.vcf.gz")


## --- map: Chr, Pos, marker --- ##
map <- data.table(
  Chr = gsub("chromosome_","Chr",VCF@fix[, "CHROM"]),
  Pos = as.integer(VCF@fix[, "POS"])
)
map[,marker:=paste(Chr,Pos,sep=":")]
## --- GTs: allelic dosage matrix (0/1/2), samples in rows, markers in columns --- ##
gt_str <- extract.gt(VCF, element = "GT")

dosage <- vapply(
  strsplit(gt_str, "[/|]"),
  function(a) {
    a <- suppressWarnings(as.integer(a))
    if (anyNA(a)) NA_integer_ else sum(a)
  },
  integer(1)
)

GTs <- matrix(dosage, nrow = nrow(gt_str), dimnames = list(map$marker, colnames(gt_str)))
GTs <- t(GTs)

saveRDS(list(GTs=GTs,map=map),"./data/parsed_vcf.rds")