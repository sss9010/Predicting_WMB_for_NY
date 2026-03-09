####################################################
# R script for analysis of all possible situations
# 1. without marker and without spatial
# 2. Without marekr but with spatial
# 3. With marker but without spatial anlaysis
# 4. with marker and with spatial analysis

###############################################
# Data preparation for the analysis
library(tidyverse)
library(sommer)
library(caret)
library(dplyr)
library(lme4)
library(rrBLUP)
library('sommer')
library(MASS)
library(dplyr)
library(tidyr)
###########################################################################################################################
# Import phenotype and genotpye data
###########################################################################################################################
library(readxl)
library(readxl)
WMBpheno <- read_excel("data/pheno/master/WMBpheno.xlsx",
                       col_types = c("text", "text", "text",
                                     "text", "text", "numeric", "numeric",
                                     "numeric", "text", "text", "numeric",
                                     "text", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "text", "text",
                                     "numeric", "numeric", "numeric",
                                     "text", "numeric", "numeric", "text",
                                     "numeric", "text"))


data<- WMBpheno %>% filter(Env %in%c("HELF24", "MCG23", "SNY22", "KET21")) %>% filter(yield_kgha > 0)

table(data$Env)
table(data$Block, data$Env)
#load("C:/Users/Siim Sepp/NY_WMB/data/wmb_multi_alldata.Rdata")
#data<- multi  %>% filter(timepoint == "1")

table(data$Row, data$Env)
###################################
#PHENO CHECK AND FILTER
#########################################
#load("data/pheno/Field_data_2020_2023.Rdata")


## **1. Importing the genotypic data**

load("data/geno/WMB_master_geno.RData")

ASRgenomics::kinship.diagnostics(A1)
#Filter the phenotypes to make sure that all of them have genomic data
indp = data$GID %in% rownames(A1) # the GIDs that have the genomic markers
table(indp)
data = data[indp,]

all(data$GID %in% rownames(A1))

#Filter the genotypes to make sure that all of them have phenotype data
indg<-  rownames(A1) %in% data$GID
table(indg)
A1<- A1[indg, indg]
dim(A1)
check<- match.kinship2pheno(K=A1, pheno.data = data, indiv = 'GID', clean = FALSE, mism = TRUE)
Ghat.bend <- G.tuneup(G= A1, bend = TRUE, eig.tol= 1e-03 )$Gb
Ginv.sparse<- G.inverse(G= Ghat.bend, sparse= TRUE)$Ginv
#ASRgenomics::kinship.heatmap(A1)
#ASRgenomics::kinship.pca(A1)


pheno<- data
str(pheno)
################################################
table(pheno$Population)
table(pheno$Env)
#pheno<- pheno[,c(1,3,4,6,7,8,9,12, 14, 21, 24, 29, 31, 36)]


## CHEC THE PHENO VARIABLES
#pheno<- pheno %>% filter(GID %in% rownames(Am))
pheno$GID = as.factor(pheno$GID)
pheno$Block = as.factor(pheno$Block)
pheno$Column = as.factor(pheno$Column)
pheno$Row = as.factor(pheno$Row)
pheno$Env = as.factor(pheno$Env)
pheno$TWT = pheno$TWT
pheno$Yld = pheno$yield_kgha
pheno$Ht = as.numeric(pheno$Ht)
pheno$HD_JD = as.numeric(pheno$HD_JD)
pheno$winter_survival = as.numeric(pheno$winter_survival)
#pheno$family = as.factor(pheno$family)
#pheno$Population = as.factor(pheno$Population)
#pheno$Check = as.factor(pheno$Check)
#pheno$SourcePLOT<- as.factor(pheno$SourcePLOT)
str(pheno)



#SPECIFY INPUT TRAIT AND ENV LEVELS
Traits = c("Yld")
Trait= 'Yld'
Env = levels(pheno$Env)

############################################################################################################################
############################################################################################################################
############################################################################################################################
# 1. The cross validation of the  without spatial and without marker
############################################################################################################################
###
# With only block row, and block column
Acc_wom_wosp = tibble()

table(pheno$Env)

#spliting the data into trainign and testing set
for(env in Env){
  SL <- subset(x = pheno, subset = Env == env) # Subseting the data for each env
  TraitN = colnames(SL[Traits])#[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))
  head(SL)
  SL$Block = as.factor(SL$Block)

  #for(Trait in TraitN){

    #Choosing the method of outlier testing for replicated and unreplicated trials
    #if(length(SL$rep)/length(levels(SL$geno)) <= 1){
      # removing outlier using boxplotstat for unreplicated trials
   #   out_ind <- which(SL[,paste(Trait)] %in% boxplot.stats(SL[,paste(Trait)])$out)

    #  if(length(out_ind) == 0){
     #   SL = SL}else{

      #    SL = SL[-out_ind,]
       # }

    #}else{
      #removing outlier for replicated trials
      #eval(parse(text = paste("outl1 <- lmer(",Trait,"~(1|GID),
       #        data= SL)")))

      #outlier = which(stats::rstudent(outl1) > 3)
      #if(length(outlier) == 0){
        #SL = SL}else{

         # SL = SL[-outlier,]
        #}
    #}
    # Creating a folder that contain 5 subset with 100 times with a total of 500
    fold5 = caret::createMultiFolds(y = unique(pheno$GID), k = 5, times = 5)


    for(i in 1:length(fold5)){
      index = fold5[[i]] # the index of the sample for training set
      #subset the phenotypic data
      train_geno = droplevels(unique(SL$GID)[index])
      train_geno_ind = which(SL$GID %in% train_geno)
      train.data <- droplevels(SL %>%
                                 filter(row_number() %in% train_geno_ind)) # subset the training set
      dim(train.data)
      test.data <- droplevels(SL %>%
                                filter(!row_number() %in% train_geno_ind)) # subset the testing set
      dim(test.data)

      #test.data[,TraitN] = NA # change the grain yield of the training set to NA value

      mod_dat = rbind(train.data, test.data) # combine the the data set for analysis

      #####################
      # make the dsisgn matrix for the blk_rwo and Blk-col
      # Random factor matrix
      #idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
      idcol <- factor(test.data$Block)
      Z.Blk.test <- model.matrix(~idcol - 1)
      rownames(Z.Blk.test) <- test.data$SourcePLOT

      eval(parse(text = paste("ans <- mmer(",Trait,"~1,
                         random=~ Block,
                         rcov=~vsr(units),
                         data= train.data)")))
      #########################
      # blockrow and blockcol effects
      befall =  as.matrix(ans$U$Block[[Trait]])
      len_b =as.numeric(levels(as.factor(test.data$Block)))
      blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set
      obs.test = test.data[,c("SourcePLOT",Trait)]
      head(test.data)
      efftest = blkeff
      r = cbind(env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs"),3))
      colnames(r)[3] = "predictability"
      Acc_wom_wosp = rbind(Acc_wom_wosp,r)

    }
  }
#RUN THE DATA INTO A DATASET

df1 = as.data.frame(Acc_wom_wosp)
head(df1)
df1$method = "Block"
df1$env<- as.factor(df1$env)
df1$predictability<- as.numeric(df1$predictability)
boxplot(df1$predictability ~ df1$env, ylim= c(-0.2,1))
title('Yld with Block')
summary(df1$predictability)



#######################################################################################################################################
# Cross validation for spatial analysis but not marker data
######################################################################################################################################
str(pheno)
pheno$GID<- as.factor(pheno$GID)
pheno$Block<- as.factor(pheno$Block)
pheno$Row<- as.numeric(pheno$Row)
pheno$Column<- as.numeric(pheno$Column)
#Traits<- 'yield_kgha'

Acc_wom_wsp = tibble()

for(env in Env){
  SL <- subset(x = pheno, subset = Env == env) # Subseting the data for each env
  TraitN = colnames(SL[Traits])#[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))
  head(SL)
  fold5 = caret::createMultiFolds(y = unique(pheno$GID), k = 5, times = 20)
  for(i in 1:length(fold5)){
    index = fold5[[i]] # the index of the sample for training set
      #subset the phenotypic data
      train_geno = droplevels(unique(SL$GID)[index])
      train_geno_ind = which(SL$GID %in% train_geno)
      train.data <- droplevels(SL %>%
                                 filter(row_number() %in% train_geno_ind)) # subset the training set
      dim(train.data)
      test.data <- droplevels(SL %>%
                                filter(!row_number() %in% train_geno_ind)) # subset the testing set
      dim(test.data)

      #test.data[,TraitN] = NA # change the grain yield of the training set to NA value

      mod_dat = rbind(train.data, test.data) # combine the the data set for analysis

      #####################
      # make the dsisgn matrix for the blk_rwo and Blk-col
      # Random factor matrix
      idcol <- factor(test.data$Block)
      #idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
      Z.Blk.test <- model.matrix(~idcol - 1)
      rownames(Z.Blk.test) <- test.data$SourcePLOT

      #add spline based row x col matrix to random effects
      eval(parse(text = paste("ans1 <- mmer(",Trait,"~1,
                             random=~ Block+ spl2Da(Row,Column),
                             rcov=~vsr(units),
                             data= train.data)")))
      befall =  as.matrix(ans1$U$Block[[Trait]])
      len_b =as.numeric(levels(test.data$Block))
      blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set

      # make a plot to observe the spatial effects found by the spl2D()
      W <- with(test.data,spl2Da(Row,Column)) # 2D spline incidence matrix
      test.data$spatial <- W$Z$`A:all`%*%ans1$U$`A:all`[[Trait]] # 2D spline BLUPs

      obs.test = test.data[,c("SourcePLOT",Trait)]
      efftest = blkeff + test.data$spatial

      r = cbind(env, Trait, predictability = round(cor(efftest[,1],obs.test[,2], use = "pairwise.complete.obs" ),3))
      colnames(r)[3] = "predictability"
      Acc_wom_wsp = rbind(Acc_wom_wsp,r)


    }
  }

#RUN THE DATA INTO THE DATASET
df2 = as.data.frame(Acc_wom_wsp)
df2$method = "Block+Spatial*"
df2$env<- as.factor(df2$env)
df2$predictability<- as.numeric(df2$predictability)
boxplot(df2$predictability ~ df2$env)
title('Yld with Block + Spatial*')
summary(df2$predictability)
dim(df2)
head(df2)


##
sp<- df2[,c(1,3)] %>% arrange(env)
dim(sp)
colnames(sp)<- c('Env', 'Spatial + Block')
#write.csv(gp, file= 'spatial.csv') EXTRACT OUT




##############################################################################################################################
# 3. Genomic selection without spatial analysis
##############################################################################################################################

#making sure the pheno data has geno data and vice versa
all(rownames(A1)%in% pheno$GID)
all(pheno$GID %in% rownames(A1))

###############################################################################################
# 4. with marker and blocrow and blkcol
############################################################################################################

Acc_wm_wosp_bkrc = tibble()
for(env in Env){
  SL <- subset(x = pheno, subset = Env == env) # Subseting the data for each env
  TraitN = colnames(SL[Traits])#[colSums(is.na(SL[Traits])) < 25] # selecting the trait
  ntt = length((TraitN))
  head(SL)
  # Creating a folder that contain 5 subset with 100 times with a total of 500
  fold5 = caret::createMultiFolds(y = unique(pheno$GID), k = 5, times = 5)
  for(i in 1:length(fold5)){
    index = fold5[[i]] # the index of the sample for training set
    #subset the phenotypic data
    train_geno = droplevels(unique(SL$GID)[index])
    train_geno_ind = which(SL$GID %in% train_geno)
    train.data <- droplevels(SL %>%
                               filter(row_number() %in% train_geno_ind)) # subset the training set
    dim(train.data)
      test.data <- droplevels(SL %>%
                                filter(!row_number() %in% train_geno_ind)) # subset the testing set
      dim(test.data)

      #test.data[,TraitN] = NA # change the grain yield of the training set to NA value
      mod_dat = rbind(train.data, test.data) # combine the the data set for analysis

      #####################
      # make the dsisgn matrix for the blk_rwo and Blk-col
      # Random factor matrix
      idcol <- factor(test.data$Block)
      #idcol <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
      Z.Blk.test <- model.matrix(~idcol - 1)
      rownames(Z.Blk.test) <- test.data$SourcePLOT

      idgtest <- unique(test.data$GID)
      #idgtest <- factor(as.character(test.data[,"GID"]), levels = unique(test.data[,"GID"]))
      Z.geno.test <- model.matrix(~idgtest - 1)
      rownames(Z.geno.test) <- levels(test.data$GID)


      eval(parse(text = paste("ans4 <- mmer(",Trait,"~1,
                                 random=~ Block + vsr(GID,Gu = A1),
                                 rcov=~vsr(units),
                                 data= train.data)")))

      genoUef = as.matrix(ans4$U$`u:GID`[[Trait]])
      genoUef= as.data.frame(genoUef)
      genoUef$GID = rownames(genoUef)
      test.data$genoeff = NA
      test.data$blockeff = NA

      befall =  as.matrix(ans4$U$Block[[Trait]])
      len_b =as.numeric(levels(test.data$Block))
      blkeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set

      test.data$blockeff =  blkeff

      for(g in as.vector(test.data$GID)){

        test.data[test.data$GID == g,"genoeff"] = genoUef[genoUef$GID == g, "V1"]
      }

      test.data$toteff = test.data[,"blockeff"]  + test.data[,"genoeff"]

      r = cbind(env, Trait, predictability = round(cor(test.data[,Trait],test.data[,"toteff"], use = "pairwise.complete.obs"),3))
      colnames(r)[3] <- c("predictability")
      Acc_wm_wosp_bkrc = rbind(Acc_wm_wosp_bkrc,r)


    }
  }

df4 = as.data.frame(Acc_wm_wosp_bkrc)
df4$method = "Block+Marker"
df4$env<- as.factor(df4$env)
df4$predictability<- as.numeric(df4$predictability)
boxplot(df4$predictability ~ df4$env,ylim= c(-0.2,1))
title('Block +Marker')


boxplot(df2$predictability ~ df2$env, ylim= c(-0.2,1))
title('Yld with Block + Spatial')


boxplot(df1$predictability ~ df1$env, ylim= c(-0.2,1))
title('Yld with Block')




###############################################################################
# 5. Spatial analysis with marker and design matrix
###############################################################################



all(pheno$GID %in% rownames(A1))
dim(A1)
dim(pheno)
library(ASRgenomics)
check<- match.kinship2pheno(K=A1, pheno.data = pheno, indiv = 'GID', clean = FALSE, mism = TRUE)




Acc_wm_wsp = tibble()

for(env in Env){
  SL <- droplevels(subset(x = pheno, subset = Env == env)) # Subseting the data for each env
  TraitN = colnames(SL[Traits])[colSums(is.na(SL[Traits])) < 25] # selecting the trait

  ntt = length((TraitN))

    fold5 = caret::createMultiFolds(y = unique(pheno$GID), k = 5, times = 5)

    for(i in 1:length(fold5)){
      index = fold5[[i]] # the index of the sample for training set
      #subset the phenotypic data
      train_geno = droplevels(unique(SL$GID)[index])
      train_geno_ind = which(SL$GID %in% train_geno)
      train.data <- droplevels(SL %>%
                                 filter(row_number() %in% train_geno_ind)) # subset the training set
      dim(train.data)
      test.data <- droplevels(SL %>%
                                filter(!row_number() %in% train_geno_ind)) # subset the testing set
      dim(test.data)

      # test.data[,TraitN] = NA # change the grain yield of the training set to NA value

      mod_dat = rbind(train.data, test.data) # combine the the data set for analysis

      #####################
      # make the dsisgn matrix for the blk_rwo and Blk-col
      # Random factor matrix
      dim(test.data)
      #idblock <- factor(as.character(test.data[,"Block"]), levels = unique(test.data[,"Block"]))
      idblock <- factor(test.data$Block)
      Z.Blk.test <- model.matrix(~idblock - 1)

      rownames(Z.Blk.test) <- test.data$SourcePLOT


      eval(parse(text = paste("ans3<- mmer(",Trait,"~1,
                                     random=~ Block + vsr(GID, Gu = A1) +
                                       spl2Da(Row,Column),
                                     rcov=~vsr(units),
                                     data= train.data)")))

      # Set the levels of test.data$Block based on unique values
      #test.data$Block <- factor(test.data$Block, levels = unique(test.data$Block))

      befall =  as.matrix(ans3$U$Block[[Trait]])
      len_b =as.numeric(levels(test.data$Block))
      blockeff = Z.Blk.test %*% befall[len_b,] # block effect for the test set

      # make a plot to observe the spatial effects found by the spl2D()
      W <- with(test.data,spl2Da(Row,Column)) # 2D spline incidence matrix
      test.data$spatial <- W$Z$`A:all`%*%ans3$U$`A:all`[[Trait]] # 2D spline BLUPs

      genoUef = as.matrix(ans3$U$`u:GID`[[Trait]])
      genoUef = as.matrix(genoUef[order(rownames(genoUef)),])
      rownames(genoUef) %in% rownames(A1)
      id = which(rownames(genoUef) %in% unique(test.data$GID))
      genoUtest = as.matrix(genoUef[id,])

      rownames(genoUtest) %in% unique(test.data$GID)
      dim(genoUtest)
      # geneff = as.matrix(ans3$U$`u:geno`[[Trait]])
      # geneff = as.matrix(geneff[order(rownames(geneff)),])
      # rownames(geneff) %in% rownames(A1)
      # id = which(rownames(geneff) %in% unique(test.data1$geno))
      # genefftest = as.matrix(geneff[id,])
      # rownames(genefftest) %in% unique(test.data1$geno)
      test.data$genoUtest = NA
      genoUtest = as.data.frame(genoUtest)
      genoUtest$GID = rownames(genoUtest)
      ################################################
      # putting the genetic effect for the indvidual plot
      #################################################
      #place the genoU blup in the table
      for(g in as.vector(test.data$GID)){

        test.data[test.data$GID == g,"genoUtest"] = genoUtest[genoUtest$GID == g, "V1"]
      }
      ############################
      # # putting genotypic effec to the data
      # for(g in as.vector(test.data$geno)){
      #
      #   test.data[test.data$geno == g,"genoeff"] = genefftest[genefftest$geno == g, "V1"]
      # }
      #

      efftest = blockeff + test.data$spatial + test.data$genoUtest
      test.data$toteff = efftest


      r = cbind(env,Trait,
                predictability = round(cor(test.data[,paste(Trait)],
                                            test.data[,"toteff"],
                                           use = "pairwise.complete.obs"),3))
      colnames(r)[3] = "predictability"
      Acc_wm_wsp = rbind(Acc_wm_wsp,r)


    }

  }


df5 = as.data.frame(Acc_wm_wsp)
df5$method = "Block+Marker+Spatial"
head(df5)
df5$env<- as.factor(df5$env)
df5$predictability <- as.numeric(df5$predictability)
boxplot(df5$predictability ~ df5$env)

df521<- df5 %>% filter(env == '21KET')
df522<- df5 %>% filter(env == '22SNY')
df523<- df5 %>% filter(env == '23MCG')
summary(df521$predictability)
summary(df522$predictability)
summary(df523$predictability)
dim(df5)
gp21<- df5 %>% filter(env == '21KET')
gp22<- df5 %>% filter(env == '22SNY')
gp23<- df5 %>% filter(env == '23MCG')

summary(gp21$predictability)
summary(gp22$predictability)
summary(gp23$predictability)

gp<- cbind(gp21$predictability, gp22$predictability, gp23$predictability)
colnames(gp)<- c('KET21', 'SNY22','MCG23')

#write.csv(df5, file= 'GPres.csv') #or the gp file


##########################################################################################################################
par(mfrow = c(2,2))
boxplot(df1$predictability ~ df1$env, ylim= c(-0.2,1))
title(' Block')

boxplot(df2$predictability ~ df2$env, ylim= c(-0.2,1))
title(' Block + Spatial')

boxplot(df4$predictability ~ df4$env,ylim= c(-0.2,1))
title('Block +Marker')

boxplot(df5$predictability ~ df5$env,ylim= c(-0.2,1))
title('Block +Marker + Spatial')

mean(df1$predictability)
mean(df2$predictability)
mean(df4$predictability)
mean(df5$predictability)


dfcomb = rbind(df1,df2,df4,df5)
str(dfcomb)
dfcomb$predictability = as.numeric(dfcomb$predictability)
table(dfcomb$method)

## stacked box plot
colnames(dfcomb)[3]<- 'Accuracy'
dfcomb$method = factor(x = dfcomb$method, levels = c("Block", "Block+Marker", "Block+Spatial*", "Block+Marker+Spatial"))
#dfcomb$method = factor(x = dfcomb$method, levels = c( "Block+Marker",  "Block+Marker+Spatial"))
#dfcomb<- dfcomb%>% filter(method %in% c( "Block+Marker+Spatial"))
dfcomb$method = as.factor(dfcomb$method)


e <- ggplot(dfcomb, aes(x = Trait, y = Accuracy)) +
  facet_grid(~env) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
e2 <- e + geom_violin(
  aes(fill = method),
  position = position_dodge(0.9)
) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#100000","#Abc111", "#BCDE2222"))


e2 + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Prediction accuracy fo spatial correction ang genomic prediction models WMB master 21-23")

library(gridExtra)
ggarrange(e2)
grid.arrange(e2)




#####################
##########################################################################################################
