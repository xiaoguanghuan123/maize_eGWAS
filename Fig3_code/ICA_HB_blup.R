
setwd('/work3/maize_eGWAS/figures/Fig3_code/')
library(lme4)
library(parallel)

################# calculate Heritability for each of the kept moduals ###########
dfA = read.csv('ICA_matrixA_sel.csv', check.names = F, row.names = 1)
cn = colnames(dfA)
idx = colnames(dfA) %in% colnames(dfA)[duplicated(colnames(dfA))]
dfA <- dfA[, idx]
colnames(dfA) = cn[idx]
length(unique(colnames(dfA))) # 219

group <- factor(colnames(dfA))
all(rownames(dfA) == names(srg) ) #TRUE

HbICA <- parallel::mclapply(seq_len(nrow(dfA)),  function(i){
  df1 <- dfA[i,] 
  df1 <- data.frame(t(df1))
  df1$gt <- group
  colnames(df1)[1] = 'gexp'
  model = lmer(gexp ~ (1|gt)  , df1)
  summary(model)
  varComp <- as.data.frame(VarCorr(model,comp='vcov'))
  sigmas <- varComp$vcov
  H2 <- sigmas[1]/sum(sigmas)
  H2
},mc.cores=8)

dfHb = data.frame(ModelName = names(srg), ModelSize = sapply(srg, length),  
                  Heritability = unlist(HbICA))
write.csv(dfHb, file = 'moduleHb_Size.csv', row.names = F, quote = F)
scatter.smooth(log2(dfHb$ModelSize), dfHb$Heritability)

#################### Calculate BLUPs for each modual for GWAS ################
dfA = read.csv('ICA_matrixA_sel.csv', check.names = F, row.names = 1)
group <- factor(colnames(dfA))

blupOutPut <- lapply(seq_len(nrow(dfA)), function(i){
  df1 <- data.frame(t(dfA[i,])) 
  df1$gt <- group
  colnames(df1)[1] = 'signal'
  model = lmer(signal ~ (1|gt)  , df1)
  varComp <- as.data.frame(VarCorr(model,comp='vcov'))
  blup <- coef(model)$gt
  blup
  #hist(blup[,1])
})

blupOutPut <- do.call(cbind,blupOutPut)
dim(blupOutPut)
colnames(blupOutPut) <- rownames(dfA) #names used for list, for a matrix, use colnames. 
blupOutPut[1:5, 1:5]

write.csv(blupOutPut, file = 'Model_Blups.csv')

for (i in 1:nrow(blupOutPut)){
  blups <- blupOutPut[i,]
  hist(blups)
}

data.frame(unique(group), levels(group))
