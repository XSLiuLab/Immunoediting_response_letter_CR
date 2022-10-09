library(dplyr)
library(ggpubr)
library(patchwork)
gppm_meta <- readRDS("~/test/gppm_meta.rds")
gppm_meta <- gppm_meta %>% filter(variant == "nonsynonymous SNV")
mt_summ <- gppm_meta %>%
  group_by(subst_type3) %>%
  summarise(neo_per=mean(HLA_aff_mean<500)) %>% ungroup()
saveRDS(mt_summ,file = "all_sim_mt.rds")


##获取突变类型
get_mut_type <- function (mutation,ref_genome){
  a <- b <- trinucleotide <- alt_reverse <- alt <- mutation_type <- ref <- reverse <- NULL
  hsgs.installed = BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = hsgs.installed)
  ref_genome = BSgenome::getBSgenome(genome = ref_genome)
  mutation$position <- as.numeric(mutation$position)
  mutation <- mutation %>% dplyr::mutate(Start = .data$position -
                                           1, End = .data$position + 1)
  extract.tbl <- data.table::setDT(mutation)
  ss = BSgenome::getSeq(x = ref_genome, names = extract.tbl[,
                                                            chromosome], start = extract.tbl[, Start], end = extract.tbl[,
                                                                                                                         End], as.character = TRUE)
  extract.tbl[, `:=`(trinucleotide, as.character(ss))]
  sub <- data.frame(x = c("T", "C", "G", "A"), y = c("A", "G",
                                                     "C", "T"), stringsAsFactors = F)
  extract.tbl[, `:=`(a = sapply(trinucleotide, function(z) {
    sub$y[which(sub$x == substr(z, 1, 1))]
  }), b = sapply(trinucleotide, function(z) {
    sub$y[which(sub$x == substr(z, 2, 2))]
  }), c = sapply(trinucleotide, function(z) {
    sub$y[which(sub$x == substr(z, 3, 3))]
  }))]
  extract.tbl[, `:=`(reverse = stringi::stri_reverse(paste(a,
                                                           b, c, sep = "")))]
  extract.tbl <- extract.tbl[, -c("a", "b", "c")]
  extract.tbl[, `:=`(alt_reverse, sapply(alt, function(z) {
    sub$y[which(sub$x == z)]
  }))]
  extract.tbl[, `:=`(mutation_type, ifelse(ref %in% c("C",
                                                      "T"), paste(substr(trinucleotide, 1, 1), "[", substr(trinucleotide,
                                                                                                           2, 2), ">", alt, "]", substr(trinucleotide, 3, 3), sep = ""),
                                           paste(substr(reverse, 1, 1), "[", substr(reverse, 2,
                                                                                    2), ">", alt_reverse, "]", substr(reverse, 3, 3),
                                                 sep = "")))]
  extract.tbl <- extract.tbl %>% select(-c("Start", "End",
                                           "trinucleotide", "reverse", "alt_reverse"))
  return(extract.tbl)
}
all_mut_tpm_not_filter <- readRDS("~/Immunoediting/data/all_mut_tpm_not_filter.rds")
colnames(all_mut_tpm_not_filter)[2] <- "chromosome"
mut_dt_wu <- NeoEnrichment::get_mutation_type(all_mut_tpm_not_filter)
TCGA_maf_sim <- readRDS("~/Immunoediting/data/TCGA_maf_sim.rds")
TCGA_maf_sim<- TCGA_maf_sim[TCGA_maf_sim$Variant_Classification=="nonsynonymous SNV",]
TCGA_maf_sim<- TCGA_maf_sim[!is.na(TCGA_maf_sim$mRNA),]
mut_dt_letter <- TCGA_maf_sim %>%
  select(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,
         Tumor_Sample_Barcode,subst_type3,mut_HLA_mean_aff,mRNA,Hugo_Symbol) %>%
  select(-End_Position)
colnames(mut_dt_letter) <- c("chromosome","position","ref","alt","sample",
                             "mut_type3","ic50","exp","gene")
mut_dt_letter <- get_mut_type(mut_dt_letter,"hg19")
saveRDS(mut_dt_wu,file = "mut_dt_wu.rds")
saveRDS(mut_dt_letter ,file = "mut_dt_letter.rds")

##相关性
all_mut_wu <- readRDS("mut_dt_wu.rds")
all_mut_wu %>%
  group_by(mutation_type)%>%
  summarise(neo_c=mean(neo=="neo")) -> all_mut_wu_mt_summ
all_mut_letter <- readRDS("~/test/mut_dt_letter.rds")
all_mut_letter %>%
  group_by(mutation_type)%>%
  summarise(neo_c=mean(ic50<500)) -> all_mut_letter_mt_summ
all_mut_summ <- left_join(
  all_mut_letter_mt_summ %>% rename(sim_neo=neo_c),
  all_mut_wu_mt_summ %>% rename(real_neo=neo_c)
)
p1 <- ggscatter(all_mut_summ, x = "sim_neo", y = "real_neo",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 0.2, label.sep = "\n", size=6),
          title = "TCGA sim VS actual mutation",xlab = "Neoantigen proportion of simulated data",
          ylab = "Neoantigen proportion of actual data"
)

mapping <- all_mut_letter %>% select(mut_type3,mutation_type) %>% distinct_all()
all_sim_mt <- readRDS("~/test/all_sim_mt.rds")
all_sim_mt <- left_join(
  all_sim_mt %>% rename(mut_type3=subst_type3),
  mapping
)
all_mut_summ <- left_join(
  all_sim_mt %>% rename(sim_neo=neo_per),
  all_mut_wu_mt_summ %>% rename(real_neo=neo_c)
)
p2 <- ggscatter(all_mut_summ, x = "sim_neo", y = "real_neo",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 0.26, label.sep = "\n",size=6),
          title = "All GPPM sim VS TCGA actual mutation",xlab = "Neoantigen proportion of simulated data",
          ylab = "Neoantigen proportion of actual data"
)

p2 + p1

############
###不表达基因突变类型的分布差异
non_exp_wu <- all_mut_wu %>%
  filter(near(tpm_exp,0)) %>%
  group_by(mutation_type) %>%
  summarise(c=n())
non_exp_letter <- all_mut_letter %>%
  filter(exp==0) %>%
  group_by(mutation_type) %>%
  summarise(c=n())
non_exp_summ <- left_join(
  non_exp_wu %>% rename(real_per=c),
  non_exp_letter %>% rename(sim_per=c)
)
non_exp_summ <- non_exp_summ %>%
  tidyr::pivot_longer(cols = c("real_per","sim_per"),names_to = "type",
                      values_to = "per")
tt <- non_exp_summ %>%
  tidyr::pivot_wider(names_from = "mutation_type",values_from = "per") %>%
  as.data.frame()
rownames(tt) <- c("Actual","Sim")
tt <- tt %>% select(-type)
tt <- as.table(as.matrix(tt))
chisq.test(tt)
##选择前10个展示
top_10_mt <- non_exp_letter %>% slice_max(order_by = c,n=10)
non_exp_summ_sub <- non_exp_summ %>% filter(mutation_type %in% top_10_mt$mutation_type)
tt <- non_exp_summ_sub %>%
  tidyr::pivot_wider(names_from = "mutation_type",values_from = "per") %>%
  as.data.frame()
rownames(tt) <- c("Actual","Sim")
tt <- tt %>% select(-type)
tt <- as.table(as.matrix(tt))
library("graphics")
mosaicplot(tt, shade = T, las=1,border = "skyblue",
           main="Mutation type distribution of non-expressed genes",
           type = "deviance",off = 3,cex.axis=1.3)

###
library(hrbrthemes)
all_mut_wu %>% as.data.frame() %>%
  mutate(index=paste(sample,chromosome,position,ref,alt,sep = ":")) %>%
  select(index,gene,tpm_exp) %>%
  distinct_all() -> all_mut_wu1
true_gene_summ <- all_mut_wu1 %>%
  group_by(gene) %>%
  summarise(mut_counts=n(),
            exp_median=median(tpm_exp,na.rm = T)) %>%
  ungroup() %>%
  mutate(exp_type = ifelse(near(exp_median,0),"Non-expressed gene","Expressed gene")) %>%
  na.omit()
p1 <- ggplot(data=true_gene_summ,aes(x=exp_type,y=log(mut_counts)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_ipsum()+
  labs(y="log(Mutation counts)",title="Actual data",x=NULL)

##sim
all_mut_letter %>% as.data.frame() %>%
  mutate(index = paste(sample,chromosome,position,ref,alt,sep = ":")) %>%
  mutate(neo=ifelse(ic50<500,"neo","not_neo")) %>%
  select(index,gene,exp,neo) %>%
  distinct_all() -> all_mut_letter1
sim_gene_summ <- all_mut_letter1 %>%
  group_by(gene) %>%
  summarise(mut_counts=n(),
            exp_median=median(exp,na.rm = T)) %>%
  ungroup() %>%
  mutate(exp_type = ifelse(exp_median==0,"Non-expressed gene","Expressed gene")) %>%
  na.omit()
p2 <- ggplot(data=sim_gene_summ,aes(x=exp_type,y=log(mut_counts)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_ipsum()+
  labs(y="log(Mutation counts)",title="Simulated data",x=NULL)

p1 + p2

