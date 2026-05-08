rm(list = ls())

pkg <- c("survival", "magrittr", "glmnet", "snow", "mlbench", "doParallel", "caret", "klaR", "UpSetR", "Boruta")

for(i in pkg){
  library(i, character.only = T)
}

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# data_process ------------------------------------------------------------

combat <- readRDS("data/figure6_input.rds")
combat <- combat[combat$batch %in% c('Bruan_RCC_pre_aPD1_combo_tpm',
                                     'Snyder_UC_pre_aPDL1',
                                     'Hugo_SKCM_pre_aPD1'), ]


# dir.create("figure8")
# 自己的基因数据
# sig_pheno <- readRDS("data/hub_genes.rds")


# # 批量循环自己的通路基因数据
# 
# # 列出当前目录下所有开头是11-hub-genes-的文件
# fdatalist=list.files('Figure4-766pathway-hubgenes-outcome','^11-hub-genes-')
# fnumlist=list.files('Figure6-111pathway-hubgenes-model-AUC-Cor-0.3-outcome-selected','^21-111pathway-')
# fnumlist1 <- stringr::str_split(fnumlist,'-',simplify = T)[,3]#只保留前面的样本序号
# fnumlist2 <- stringr::str_split(fnumlist,'-',simplify = T)[,4]#只保留前面的样本名称
# fdatalist <- fdatalist[c(as.numeric(fnumlist1)   )]
# dir.create("Figure8-33pathway-hubgenes-6modelalgorithm-last3-outcome")
# dir.create("Figure8-33pathway-hubgenes-6modelalgorithm-outcome")
# 
# 
# # for (o in c(1:length(fdatalist))) {
# for (o in c(25:length(fdatalist))) {
#   # o <- 15
#   print(paste0( o  ," " ,fdatalist[o] ," is start!"))
#   print(Sys.time() )
#   hubgenedata <- data.table::fread(paste0("Figure4-766pathway-hubgenes-outcome/",fdatalist[o]),header = T,data.table = F)
#   sig_pheno  <- hubgenedata[,2]
#   #运行加载好的代码，得到结果
# 
# pre_var <- intersect(colnames(combat), sig_pheno)
# 
# df <- subset(combat, select = c('response', pre_var)) %>% na.omit()
# 
# df$response <- factor(df$response, levels = c('R','NR'))
# 
# # logistic-lasso ----------------------------------------------------------
# 
# source("code/figure8_function.R")
# 
# cutoff_line <- 0.5
# 
# # logistic筛选基因
# Logoutput <- NULL
# for(i in 2:ncol(df)){
#   # i <- 2
#   g <- colnames(df)[i]
#   tryCatch(
#   mod1 <- glm(response~get(colnames(df)[i]), family = binomial(link = 'logit'),data = df)
#   , error = function(err) {
#     print(paste(o, " 异常"))
#   })
#   fit <- summary(mod1)
#   Logoutput=rbind(Logoutput,data.frame(gene=g,
#                                        OR=as.numeric(fit$coefficients[,"Estimate"])[2],
#                                        z=as.numeric(fit$coefficients[,"z value"])[2],
#                                        pvalue=as.numeric(fit$coefficients[,"Pr(>|z|)"])[2],stringsAsFactors = F))
# }
# 
# 
# (log.res <- Logoutput[which(Logoutput$pvalue < 0.5),"gene"])
# 
# df <- subset(df,select = c('response',log.res))
# 
# # lassso
# tryCatch(
# diagnosis_lasso <- lasso_iter(log.res, group = "response", df, lambda_choose = 'lambda.min', nfold = 5, iter.times=1000)
# , error = function(err) {
#   print(paste(o, " 异常"))
# })
# 
# 
# 
# lasso_hubgenes <- diagnosis_lasso %>% dplyr::filter(len!=0) %>% dplyr::group_by(marker) %>% 
#   dplyr::count(genes) %>% tail(1) %>% as.data.frame() %>% dplyr::select(genes) %>%
#   stringr::str_split(., "\\|") %>% '[['(1) %>% trimws()
# 
# 
# 
# # machine learn -----------------------------------------------------------
# cutoff_line <- 0.5
# 
# tryCatch(
# ls_res <- machine_res(df)
# 
# , error = function(err) {
#   print(paste(o, " 异常"))
# })
# 
# 
# # Upset -------------------------------------------------------------------
# 
# listInput <- ls_res$hubgenes
# listInput[["lasso"]] <- lasso_hubgenes
# # hub genes ---------------------------------------------------------------
# 
# hub_genes <- unlist(listInput)
# # hub_genes2 <- c(unlist(listInput) ,c(lasso_hubgenes)  )
# hub_genes2 <- which(table(hub_genes) >= 3) %>% names
# hub_genes <- as.data.frame(hub_genes)
# colnames(hub_genes) <- fnumlist2[o]
# 
# write.csv(hub_genes,paste0("Figure8-33pathway-hubgenes-6modelalgorithm-outcome/31-33pathway-",o,"-",
#                            fnumlist2[o],"-6modelalgorithm-outcome.csv"))
# 
# hub_genes2 <- as.data.frame(hub_genes2)
# colnames(hub_genes2) <- fnumlist2[o]
# 
# write.csv(hub_genes2,paste0("Figure8-33pathway-hubgenes-6modelalgorithm-last3-outcome/31-33pathway-",o,"-",
#                            fnumlist2[o],"-6modelalgorithm-last3-",nrow(hub_genes2),"-outcome.csv"))
# 
# print(paste0( o  ," " ,fdatalist[o] ," is end!"))
# print(Sys.time() )
# 
# }
# 
# 
# 
# listInput[["lasso"]] <- lasso_hubgenes
# 
# pdf(file = 'figure8/figG.pdf',height = 5,width = 6,onefile = F)
# upset(fromList(listInput), order.by = "freq",nsets = length(listInput))
# dev.off()





# sig_pheno <- readRDS("data/hub_genes.rds")

pheno_genes <- readRDS("figure4/1-E2F1_Q3_01-pathway-gene-Cor-0.3-pvalue-0.05-selected.rds")
sig_pheno<- pheno_genes$E2F1_Q3_01


pre_var <- intersect(colnames(combat), sig_pheno)

df <- subset(combat, select = c('response', pre_var)) %>% na.omit()

df$response <- factor(df$response, levels = c('R','NR'))

# logistic-lasso ----------------------------------------------------------

source("code/figure8_function.R")

cutoff_line <- 0.5

# logistic筛选基因
Logoutput <- NULL
for(i in 2:ncol(df)){
  # i <- 2
  g <- colnames(df)[i]
  mod1 <- glm(response~get(colnames(df)[i]), family = binomial(link = 'logit'),data = df)
  fit <- summary(mod1)
  Logoutput=rbind(Logoutput,data.frame(gene=g,
                                       OR=as.numeric(fit$coefficients[,"Estimate"])[2],
                                       z=as.numeric(fit$coefficients[,"z value"])[2],
                                       pvalue=as.numeric(fit$coefficients[,"Pr(>|z|)"])[2],stringsAsFactors = F))
}


(log.res <- Logoutput[which(Logoutput$pvalue < 0.5),"gene"])

df <- subset(df,select = c('response',log.res))

# lassso

diagnosis_lasso <- lasso_iter(log.res, group = "response", df, lambda_choose = 'lambda.min', nfold = 5, iter.times=1000)

lasso_hubgenes <- diagnosis_lasso %>% dplyr::filter(len!=0) %>% dplyr::group_by(marker) %>% 
  dplyr::count(genes) %>% tail(1) %>% as.data.frame() %>% dplyr::select(genes) %>%
  stringr::str_split(., "\\|") %>% '[['(1) %>% trimws()



# machine learn -----------------------------------------------------------
cutoff_line <- 0.5
ls_res <- machine_res(df)


dir.create("figure8")
# Upset -------------------------------------------------------------------

listInput <- ls_res$hubgenes
listInput[["lasso"]] <- lasso_hubgenes

pdf(file = 'figure8/figure8-figG-lasso.pdf',height = 5,width = 6,onefile = F)
upset(fromList(listInput), order.by = "freq",nsets = length(listInput))
dev.off()





# plot --------------------------------------------------------------------

#bar_plot-lasso
counts <- diagnosis_lasso %>% dplyr::count(marker)
counts <- counts[-1,]
diagnosis_lasso_plot <- subset(diagnosis_lasso, marker != "0 genes A")

E <- ggplot()+
  geom_bar(data = diagnosis_lasso_plot, mapping = aes(marker), fill = "#00CDCD")+
  ylab("Frequency")+
  xlab("")+
  ggtitle("Frequency of models")+
  geom_text(data = counts, aes(label = n, x = marker, y = n), vjust = -0.5)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black",angle = 90),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))
ggsave("figure8/figure8-figE-bar_plot-lasso.pdf", plot = E, height = 6, width = 6)

# lqv
C <- ggplot(ls_res$res$lqv_res) + geom_hline(aes(yintercept= cutoff_line),col='red',lty=3) + 
  ggtitle("LQV")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figure8/figure8-figC-LQV.pdf", plot = C, height = 6, width = 4)

# Bagged trees
pdf(file = 'figure8/figure8-figA-Bagged trees.pdf',height = 6,width = 4)
plot(ls_res$res$treebag_res, type=c("o"), main = "Bagged trees")
dev.off()

# Boruta
pdf(file = 'figure8/figure8-figD-Boruta.pdf',height = 6,width = 8)
plot(ls_res$res$boruta_res, las=2, xlab = "", cex.axis = 0.8, main = "Boruta")
dev.off()

# Bayesian
pdf(file = 'figure8/figure8-figB-Bayesian.pdf',height = 6,width = 4)
plot(ls_res$res$bayesian_res, type=c("o"), main = "Bayesian")
dev.off()

# Random Forest
pdf(file = 'figure8/figure8-figF-Random Forest.pdf',height = 4,width = 6)
plot(ls_res$res$rfe_res, type=c("o"), main = "Random Forest")
dev.off()


# hub genes ---------------------------------------------------------------

hub_genes <- unlist(listInput)

hub_genes <- which(table(hub_genes) >= 3) %>% names

saveRDS(hub_genes, file = 'data/hub_genes_final.rds')

