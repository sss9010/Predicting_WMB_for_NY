library(readxl)
library(caret)
library(glmnet)
library(mlbench)
library(psych)

load(file= 'data/pheno/master/flightdataWMB21_23.Rdata')
head(spctr)
colnames(spctr)

table(spct$Env)

colnames(pheno)
pi<- all_pheno[,c(1,6,7,8)]

spctry<- merge(pi, spctr, by= 'SourcePLOT')
table(spctry$Env)
# Assuming spctr is your data frame with columns 'red', 'green', 'blue', 'nir', and 'red_edge'

# Calculate vegetation indices directly in the 'spctr' data frame
#spctr$BCC <- spctr$blue / (spctr$red + spctr$green + spctr$blue)
#spctr$BGI <- spctr$blue / spctr$green
#spctr$BI <- sqrt((spctr$red^2 + spctr$green^2 + spctr$blue^2) / 3)
#spctr$CIVE <- 0.441 * spctr$red - 0.811 * spctr$green + 0.385 * spctr$blue + 18.78745
#spctr$COM1 <- spctr$EXG + spctr$CIVE + spctr$EXGR + spctr$VEG
#spctr$COM2 <- 0.36 * spctr$EXG + 0.47 * spctr$CIVE + 0.17 * spctr$VEG
#spctr$EBI <- (spctr$blue - spctr$green) / (spctr$blue - spctr$red)
#spctr$EGI <- (spctr$green - spctr$red) / (spctr$red - spctr$blue)
#spctr$ERI <- (spctr$red - spctr$green) / (spctr$red - spctr$blue)
#spctr$EXG <- 2 * spctr$green - spctr$red - spctr$blue
#spctr$EXG2 <- (2 * spctr$green - spctr$red - spctr$blue) / (spctr$green + spctr$red + spctr$blue)
#spctr$EXGR <- 3 * spctr$green - 2.4 * spctr$red - spctr$blue
#spctr$EXR <- 1.4 * spctr$red - spctr$green
#spctr$GB <- spctr$green - spctr$blue
#spctr$GR <- spctr$green - spctr$red
#spctr$GB_SR <- spctr$green / spctr$blue
#spctr$GR_SR <- spctr$green / spctr$red
#spctr$GCC <- spctr$green / (spctr$red + spctr$green + spctr$blue)
#spctr$GLI <- (2 * spctr$green - spctr$red - spctr$blue) / (2 * spctr$green + spctr$red + spctr$blue)
#spctr$MEXG <- 1.262 * spctr$green - 0.884 * spctr$red - 0.311 * spctr$blue
#spctr$MGVRI <- (spctr$green^2 - spctr$red^2) / (spctr$green^2 + spctr$red^2)
#spctr$NDI <- spctr$nir - spctr$red_edge
#spctr$ENDVI <- (spctr$nir + spctr$green - 2 * spctr$blue) / (spctr$nir + spctr$green + 2 * spctr$blue)
#spctr$EVI <- 2.5 * (spctr$nir - spctr$red) / (spctr$nir + 6 * spctr$red - 7.5 * spctr$blue + 1)
#spctr$GDVI <- spctr$nir - spctr$green
#spctr$GIPVI <- spctr$nir / (spctr$nir + spctr$green)
#spctr$GNDVI <- (spctr$nir - spctr$green) / (spctr$nir + spctr$green)
#spctr$GOSAVI <- (1 + 0.16) * (spctr$nir - spctr$green) / (spctr$nir + spctr$green + 0.16)
#spctr$GRDVI <- (spctr$nir - spctr$green) / sqrt(spctr$nir + spctr$green)
#spctr$GRVI <- spctr$nir / spctr$green
#spctr$GSAVI <- 1.5 * ((spctr$nir - spctr$green) / (spctr$nir + spctr$green + 0.5))
#spctr$GWDRVI <- (0.12 * spctr$nir - spctr$green) / (0.12 * spctr$nir + spctr$green)
#spctr$MDD <- (spctr$nir - spctr$red_edge) - (spctr$red_edge - spctr$green)
#spctr$MGSAVI <- 0.5 * (2 * spctr$nir + 1 - sqrt((2 * spctr$nir + 1)^2 - 8 * (spctr$nir - spctr$green)))
#spctr$MNDI <- (spctr$nir - spctr$red_edge) / (spctr$nir - spctr$green)
#spctr$MNDRE <- (spctr$nir - (spctr$red_edge - 2 * spctr$green)) / (spctr$nir + (spctr$red_edge - 2 * spctr$green))
##spctr$MRESAVI <- 0.5 * (2 * spctr$nir + 1 - sqrt((2 * spctr$nir + 1)^2 - 8 * (spctr$nir - spctr$red_edge)))
#spctr$MRETVI <- 1.2 * (1.2 * (spctr$nir - spctr$green) - 2.5 * (spctr$red_edge - spctr$green))
#spctr$MSR <- ((spctr$nir / spctr$red - 1)) / sqrt((spctr$nir / spctr$red + 1))
#spctr$MSR_G <- ((spctr$nir / spctr$green - 1)) / sqrt((spctr$nir / spctr$green + 1))
#spctr$MSR_RE <- ((spctr$nir / spctr$red_edge - 1)) / sqrt((spctr$nir / spctr$red_edge + 1))
#spctr$MTCARI <- 3 * ((spctr$nir - spctr$red_edge) - 0.2 * (spctr$nir - spctr$green) * (spctr$nir / spctr$red_edge))
spctr$NDRE <- (spctr$nir - spctr$red_edge) / (spctr$nir + spctr$red_edge)
spctr$NDVI <- (spctr$nir - spctr$red) / (spctr$nir + spctr$red)
#spctr$NNIR <- spctr$nir / (spctr$nir + spctr$red_edge + spctr$green)
#spctr$NREI <- spctr$red_edge / (spctr$nir + spctr$red_edge + spctr$green)
#spctr$NGI <- spctr$green / (spctr$nir + spctr$red_edge + spctr$green)
#spctr$OSAVI <- (spctr$nir - spctr$red) / (spctr$nir + spctr$red + 0.16)
#spctr$PSRI <- (spctr$red - spctr$green) / spctr$red_edge
#spctr$REGDVI <- spctr$red_edge - spctr$green
#spctr$REGNDVI <- (spctr$red_edge - spctr$green) / (spctr$red_edge + spctr$green)
#spctr$REGRVI <- spctr$red_edge / spctr$green
#spctr$REOSAVI <- (1 + 0.16) * (spctr$nir - spctr$red_edge) / (spctr$nir + spctr$red_edge + 0.16)
##spctr$RERDVI <- (spctr$nir - spctr$red_edge) / sqrt(spctr$nir + spctr$red_edge)
#spctr$RESAVI <- 1.5 * ((spctr$nir - spctr$red_edge) / (spctr$nir + spctr$red_edge + 0.5))
#spctr$RETVI <- 0.5 * (120 * (spctr$nir - spctr$green) - 200 * (spctr$red_edge - spctr$green))
#spctr$REWDRVI <- (0.12 * spctr$nir - spctr$red_edge) / (0.12 * spctr$nir + spctr$red_edge)
#spctr$RVI <- spctr$nir / spctr$red
#spctr$SAVI <- 1.5 * ((spctr$nir - spctr$red) / (spctr$nir + spctr$red + 0.5))
#spctr$TVI <- 0.5 * (120 * (spctr$nir - spctr$green) - 200 * (spctr$red - spctr$green))
#spctr$VIopt1 <- 100 * (log(spctr$nir) - log(spctr$red_edge))
#spctr$TNDVI <- sqrt((spctr$nir - spctr$red) / (spctr$nir + spctr$red) + 0.5)
#spctr$MNLI <- 1.5 * ((spctr$nir^2 - spctr$red) / (spctr$nir^2 + spctr$red + 0.5))
#spctr$RESR <- spctr$red_edge / spctr$red
#spctr$RENDVI <- (spctr$red_edge - spctr$red) / (spctr$red_edge + spctr$red)
#spctr$NNIR2 <- spctr$nir / (spctr$nir + spctr$red_edge + spctr$red)
#spctr$NREI2 <- spctr$red_edge / (spctr$nir + spctr$red_edge + spctr$red)
#spctr$NRI <- spctr$red / (spctr$nir + spctr$red_edge + spctr$red)
ASRgenomics::kinship.heatmap(K=cor(spctr[,c(-9:-1)], use= 'pairwise.complete.obs'))
str(all_pheno)
#######################################################################################################################
all_pheno<- merge(all_pheno,spctr , by= 'SourcePLOT')
spctry<- all_pheno
#MCG23<- merge(MCG23, NDVI23, by= 'SourcePLOT')
#SNY22<- merge(SNY22, NDVI22, by= 'SourcePLOT')
#KET21<- merge(KET21, NDVI21, by= 'SourcePLOT')
#ndvi, savi, evi, evi, ndwi, !!!chi, tgi, ndre
spctr<- spctr %>% na.omit()

str(spctr)

spc21<- spctr %>% filter(Env == 'KET21')

spc22<- spctr %>% filter(Env == 'SNY22')
colnames(spc22)
spc23<- spctr %>% filter(Env == 'MCG23')



write.csv(spctr, 'NDVIRE.csv')
##############################################################################################################################################################
#Pivoting the data from TP columns  longer. 
head(spc21ndvi)
spc21ndvi<- spc21 %>% pivot_wider(
  id_cols= c(Env, yield_kgha, SourcePLOT),
  names_from = JD, 
  values_from = NDVI)

spc22ndvi<- spc22 %>% pivot_wider(
  id_cols= c(Env, yield_kgha, SourcePLOT),
  names_from = JD, 
  values_from = NDVI)

spc23ndvi<- spc23 %>% pivot_wider(
  id_cols= c(Env, yield_kgha, SourcePLOT),
  names_from = JD, 
  values_from = NDVI)



round(cor(spc[,c(2,5:8)], use= 'complete'),2)
ASRgenomics::kinship.heatmap(K=cor(spc[,c(2,5:8)], use= 'pairwise.complete'))
#ndvi, savi, evi, evi, ndwi, !!!chi, tgi, ndre



ggplot(NDlong, aes(x=Timepoint, y=NDVI_value)) +
  geom_boxplot(aes(color= Timepoint))+
  scale_color_manual(values = color_palette) +
  theme(legend.position = "none")







##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#ML based phenomic preiciton models for 3 years of WMB and 2 indices. 
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

############### Required Libraries ########################



i<- spc21ndre[,c(-1,-3)] %>% na.omit()
i<- spc22ndre[,c(-1, -3)] %>% na.omit()
i<- spc23ndre[,c(-1, -3)] %>% na.omit()
i<- spc21ndvi[,c(-1, -3)] %>% na.omit()
i<- spc22ndvi[,c(-1, -3)] %>% na.omit()
i<- spc23ndvi[,c(-1, -3)] %>% na.omit()

traits =1
cycles = 100
accuracy_lm = matrix(nrow = cycles, ncol = traits)
accuracy_ridge = matrix(nrow = cycles, ncol = traits)
accuracy_lasso = matrix(nrow = cycles, ncol = traits)
accuracy_en = matrix(nrow = cycles, ncol = traits)
accuracy_rf = matrix(nrow = cycles, ncol = traits)
#accuracy_nnet = matrix(nrow = cycles, ncol = traits)

lmI<- data.frame()
ridgeI<- data.frame()
lassoI<- data.frame()
enI<- data.frame()
rfI<- data.frame()

colnames(i)

#save(importanceNDVI, file = 'NDVIimp.csv')


for (r in 1:cycles) 
{
  ind <- createDataPartition(i$yield_kgha, p=0.8, list=FALSE)
  train <- i[ind,]
  test <- i[-ind,]
  
  
  ################### K-Fold cross validation ################
  
  custom <- trainControl(method= "repeatedcv",
                         number=5,
                         repeats = 1,
                         verboseIter = T)
  
  ############################################################
  ###################### Linear Regression ###################
  ############################################################
 table(is.na(train$yield_kgha)) 
  
  lm <- train(yield_kgha ~.,
              train,
              method= 'lm',
              trControl= custom)
  
  summary(lm)
  
  #plot(lm$finalModel)
  
  lmI_iteration <- data.frame( Importance = varImp(lm, scale = TRUE)$importance)
  
  # Combine the results with the previous iterations
  if (r == 1) {
    lmI <- lmI_iteration
    
  } else {
    lmI <- cbind(lmI, lmI_iteration)
    
  }




  predicted_test_lm <- predict(lm,test)
  accuracy_lm[r,1] <- cor(predicted_test_lm , test$yield_kgha, use="complete")
  #cor(predicted_test_lm , test$yield_kgha, use="complete")
  
  
  ############################################################
  ###################### Ridge Regression ####################
  ############################################################
  
  ridge <- train(yield_kgha~.,
                 train,
                 method = 'glmnet',
                 tuneGrid= expand.grid(alpha=0,
                                       lambda = seq(0.0001,1, length=5)),
                 trControl= custom)
  
  
  
  #plot(ridge$finalModel, xvar = "lambda", label = T)
  #plot(ridge$finalModel, xvar = "dev", label = T)
  
  #plot(varImp(ridge,scale = F))
  plot(varImp(ridge,scale = T))
  
  
  predicted_test_ridge <- predict(ridge,test)
  accuracy_ridge[r,1] <- cor(predicted_test_ridge , test$yield_kgha, use="complete")
  cor(predicted_test_ridge , test$yield_kgha, use="complete")
  
  ridgeI_iteration <- data.frame( Importance = varImp(ridge, scale = TRUE)$importance)
  
  # Combine the results with the previous iterations
  if (r == 1) {
    ridgeI <- ridgeI_iteration
    
  } else {
    ridgeI <- cbind(ridgeI, ridgeI_iteration)
    
  }
  
  
  ############################################################
  ###################### Lasso Regression ####################
  ############################################################
  
  lasso <- train(yield_kgha~.,
                 train,
                 method = 'glmnet',
                 tuneGrid= expand.grid(alpha=1,
                                       lambda = seq(0.0001,1, length=5)),
                 trControl= custom)
  #plot(lasso)
  #plot(lasso$finalModel, xvar = "lambda", label = T)
  #plot(lasso$finalModel, xvar = "dev", label = T)
  
  plot(varImp(lasso,scale = F))
  
 
  #plot(varImp(lasso,scale = T))
  
  
  predicted_test_lasso <- predict(lasso,test)
  accuracy_lasso[r,1] <- cor(predicted_test_lasso , test$yield_kgha, use="complete")
  cor(predicted_test_lasso , test$yield_kgha, use="complete")
  
  
  
  
  lassoI_iteration <- data.frame( Importance = varImp(lasso, scale = TRUE)$importance)
  
  # Combine the results with the previous iterations
  if (r == 1) {
    lassoI <- lassoI_iteration
    
  } else {
    lassoI <- cbind(lassoI, lassoI_iteration)
    
  }
  
  
  ############################################################
  ################### Elastic Net Regression #################
  ############################################################
  
  set.seed(1234)
  
  en <- train(yield_kgha~.,
              train,
              method = 'glmnet',
              tuneGrid= expand.grid(alpha=seq(0.1, length=10),
                                    lambda = seq(0.0001,1, length=5)),
              trControl= custom)
  
  plot(en)
  #plot(en$finalModel, xvar = "lambda", label = T)
  #plot(en$finalModel, xvar = "dev", label = T)
  
  plot(varImp(en,scale = F))
 
  #plot(varImp(en,scale = T))
  
  
  predicted_test_en <- predict(en,test)
  accuracy_en[r,1] <- cor(predicted_test_en , test$yield_kgha, use="complete")
  cor(predicted_test_en , test$yield_kgha, use="complete")
  
  
  
  enI_iteration <- data.frame( Importance = varImp(en, scale = TRUE)$importance)
  
  # Combine the results with the previous iterations
  if (r == 1) {
    enI <- enI_iteration
    
  } else {
    enI <- cbind(enI, enI_iteration)
    
  }
  

  
  ############################################################
  ################### Random f#################
  ############################################################
  rf <- train(yield_kgha~.,
                    train,
                    method = 'rf',
                    trControl= custom)
  #plot(rf)
  #plot(rfl$finalModel, xvar = "lambda", label = T)
  #plot(rf$finalModel, xvar = "dev", label = T)
  
  plot(varImp(rf,scale = F))
  
  #plot(varImp(rf_model,scale = T))
  
  
  predicted_test_rf <- predict(rf,test)
  accuracy_rf[r,1] <- cor(predicted_test_rf , test$yield_kgha, use="complete")
  cor(predicted_test_rf , test$yield_kgha, use="complete")
  
  
  rfI_iteration <- data.frame( Importance = varImp(rf, scale = TRUE)$importance)
  
  # Combine the results with the previous iterations
  if (r == 1) {
    rfI <- rfI_iteration
    
  } else {
    rfI <- cbind(rfI, rfI_iteration)
    
  }
  
  
  
}
  ####################################
 
summary(spc21$JD)

lmI$mod<- 'LM'
ridgeI$mod<- 'RR'
lassoI$mod<- 'LASSO'
enI$mod<- 'EN'
rfI$mod<- 'RF'
head(importance21)

importance21<- rbind(lmI, ridgeI, lassoI, enI, rfI)
importance21$Env<- 'KET21' 
importance22<- rbind(lmI, ridgeI, lassoI, enI, rfI)
importance22$Env<- 'SNY22' 
#importance23<- rbind(lmI, ridgeI, lassoI, enI, rfI)
importance23$Env<- 'MCG23' 
str(importance21)
print(importance22)
head(importance22)
importanceNDVI<- rbind(importance21, importance22, importance23)
#importanceNDRE<- rbind(importance21, importance22, importance23)
################### Prediction accuraies of the models #####################
write.csv(importance21, file= 'NDREimp21.csv')
#importanceNDRE<- importanceNDRE %>% arrange(Env)
importanceNDVI<- importanceNDVI %>% arrange(Env)
importanceNDRE$Env
importanceNDVI$Env

head(importanceNDRE)
colnames(importanceNDRE)<- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'mod', 'Env')
it<- importanceNDRE %>% filter(Env == 'KET21')
ii<- as.data.frame(t(importanceNDRE))
ii21<- ii[,1:12]
dim(importance23)



####################################################################################################33

dim(ket21NDRE)

summary(accuracy_lm)
summary(accuracy_en)
summary(accuracy_ridge)#barely better
summary(accuracy_lasso)
summary(accuracy_rf)
accuracy_lm
#ket21NDRE<- cbind(accuracy_lm,accuracy_en, accuracy_ridge,accuracy_lasso, accuracy_rf)
#sny22NDRE<- cbind(accuracy_lm,accuracy_en, accuracy_ridge,accuracy_lasso, accuracy_rf)
#mcg23NDRE<- cbind(accuracy_lm,accuracy_en, accuracy_ridge,accuracy_lasso, accuracy_rf)

ket21NDVI<- cbind(accuracy_lm,accuracy_en, accuracy_ridge,accuracy_lasso, accuracy_rf)
sny22NDVI<- cbind(accuracy_lm,accuracy_en, accuracy_ridge,accuracy_lasso, accuracy_rf)
#mcg23NDVI<- cbind(accuracy_lm,accuracy_en, accuracy_ridge,accuracy_lasso, accuracy_rf)
head(mcg23NDVI)
mcg23NDVI
ket21NDRE<- ket21NDRE[-1]

sny22NDRE$ENV<- 'SNY22'
mcg23NDRE$ENV<- 'MCG23'


dfcmb<- rbind(ket21rest, sny22rest, mcg23rest)
dfcmb$ENV<- as.factor(dfcmb$ENV)
boxplot(mcg23resM, ylim= c(0.5,1))



colnames(ket21NDRE)<- c('LM','EN' , 'RR', 'LASSO', 'RF')
colnames(sny22NDRE)<- c('LM','EN' , 'RR', 'LASSO', 'RF')
colnames(mcg23NDRE)<- c('LM','EN' , 'RR', 'LASSO', 'RF')

colnames(ket21NDVI)<- c('LM','EN' , 'RR', 'LASSO', 'RF')
colnames(sny22NDVI)<- c('LM','EN' , 'RR', 'LASSO', 'RF')
colnames(mcg23NDVI)<- c('LM','EN' , 'RR', 'LASSO', 'RF')
par(mfrow= c(1,3))


#NDVI<- rbind(ket21NDVI, sny22NDVI, mcg23NDVI)


NDRE<- rbind(ket21NDRE, sny22NDRE, mcg23NDRE)
boxplot(ket21NDRE)
boxplot(sny22NDRE)
boxplot(mcg23NDRE)

boxplot(ket21NDVI)
boxplot(sny22NDVI)
boxplot(mcg23NDVI)




KET21NDRE<- as.data.frame(ket21NDRE) %>% mutate(Env ='KET21')
SNY22NDRE<- as.data.frame(sny22NDRE)%>% mutate(Env ='SNY22')
MCG23NDRE<- as.data.frame(mcg23NDRE)%>% mutate(Env ='MCG23')
KET21NDVI<- as.data.frame(ket21NDVI)%>% mutate(Env ='KET21')
SNY22NDVI<- as.data.frame(sny22NDVI)%>% mutate(Env ='SNY22')
MCG23NDVI<- as.data.frame(mcg23NDVI)%>% mutate(Env ='MCG23')




##################################################################
#pushing all data into long format and the gneomic prediction data
#pivoting together different indices, models, and environments. . 
NDRE21<- KET21NDRE[,-5] %>%  pivot_longer(
  cols = c("LM" ,   "EN" ,   "RR"  ,  "LASSO"  ),    # Columns to pivot
  names_to = "Model",        # New column for time points
  values_to = "Accuracy"       # New column for NDVI values
)


colnames(KET21NDRE[,-5])

colnames(MCG23NDRE)


aov<- aov(Accuracy ~  Model, data= NDRE21)
NDRE21$Model<- as.factor(NDRE21$Model)
TukeyHSD(aov)

NDRE<- rbind(KET21NDRE, SNY22NDRE, MCG23NDRE)
NDVI<- rbind(KET21NDVI, SNY22NDVI, MCG23NDVI)
table(NDRE$Env)


save(NDRE, NDVI, file= 'ndrendvi.Rdata')
load('ndrendvi.Rdata')

save(gp, sp, file= 'genpred+spat.Rdata')

gspNDVI<- cbind(gp, sp[,2], NDVI)
gspNDRE<- cbind(gp, sp[,2], NDRE)




gspNDVI$I<- 'NDVI'
gspNDRE$I<- 'NDRE'

################################################
colnames(NDVIM)
table(NDVIM$Model)


dim(gspNDRE)
NDREM<- gspNDRE[, -9] %>%  pivot_longer(
  cols = c('Spatial + Block' ,'GP' , "LM" ,   "EN" ,   "RR"  ,  "LASSO", 'RF' ),    # Columns to pivot
  names_to = "Model",        # New column for time points
  values_to = "Accuracy"       # New column for NDVI values
)

NDVIM<- gspNDVI[,-9] %>%  pivot_longer(
  cols = c('Spatial + Block' ,'GP' , "LM" ,   "EN" ,   "RR"  ,  "LASSO", 'RF' ),    # Columns to pivot
  names_to = "Model",        # New column for time points
  values_to = "Accuracy"       # New column for NDVI values
)

NDREM$Env<- as.factor(NDREM$Env)
NDREM$model<- as.factor(NDREM$Model)
NDREM$Accuracy<- as.numeric(NDREM$Accuracy)
NDVIM$Env<- as.factor(NDVIM$Env)
NDVIM$model<- as.factor(NDVIM$Model)
NDVIM$Accuracy<- as.numeric(NDVIM$Accuracy)



table(NDVIM$Model)

###############################################
########################final dataframe. 
ii<- aov( NDVIM$Accuracy~NDVIM$Env*NDVIM$Model)
TukeyHSD(ii)
ii<- aov( NDREM$Accuracy~NDREM$Env +NDREM$Model +NDREM$Env:NDREM$Model)
TukeyHSD(ii)


#####
VIRE<- rbind(NDVIM, NDREM)
str(VIRE)
VIRE$Model<- as.factor(VIRE$Model)
table(VIRR$Model)
#save(VIRE, file= 'VIRE.Rdata')
summary(VIRE$Model)
VIRE$I<- as.factor(VIRE$I)



VIRR<- VIRE %>% filter(I == 'NDVI')
VIRRE<- VIRE %>% filter(I == 'NDRE')

#ggplot(VIRR, aes(x= Env, y= Accuracy, z= Model))+
table(VIRE$Model)
VIRR$Model <- factor(VIRR$Model, levels = c('Spatial + Block', "GP", "LM", "RR", "LASSO", "EN", 'RF'))

VIRRE$Model <- factor(VIRRE$Model, levels = c('Spatial + Block',"GP", "LM", "RR", "LASSO", "EN", 'RF'))

# Create the violin plot with column labels

g1<- ggplot(VIRR, aes(x = Model, y = Accuracy, fill= Model)) +
  facet_grid(~Env)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_violin() +
  geom_crossbar(stat = "summary", fun = "mean", width = 0.5, aes(ymax = ..y.., ymin = ..y..)) +
  labs(fill = "Model")+ 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Prediction accuracy of Spatial, GP and NDRE phenomic prediction models")
  

grid.arrange(nrow=1, ncol=2, g1, g2)





############################################################################################################
#different years have different flight TPs, need to switch Tp to a random variable? 

#if all NDVI tp-s then  then 0.67 mean
# if onthe best r and prediction accuracy with almost all included variables. 

spctrndvi<- spctry %>% pivot_wider(
  id_cols= c(Env, yield_kgha.x, SourcePLOT, Block, Row, Column),
  names_from = JD, 
  values_from = NDVI)

colnames(spctry)
spct<- spctry[,c(1,2,3,4,5,11,12,13,14)]
colnames(spct)
table(spct$Env)
spct21<- spct %>% filter(Env == 'KET21') %>% na.omit()
str(spct21)

spct21$JD <- as.factor(spct21$JD)
spct21$Block<- as.factor(spct21$Block)
spct21$Row<- as.numeric(spct21$Row)
spct21$Column<- as.numeric(spct21$Column)

############################################################################################################################################################
################wEIRD COMBINED MODEL
###########################################################################################################################################################
mixmph22<- mmer(yield_kgha~JD,
              random=~ Block + JD:NDVI +
                spl2Da(Row,Column),
              rcov=~vsr(units), 
              data = spct22)


mm21<-plot(mixmph)
ipl21<- summary(mixmph21)
ipl$varcomp[2,2]/sum(ipl$varcomp[,2]) #0.92 explained by JD:NDVI
ipl$varcomp
ipl$betas
ipl$logo
summary(lm)

###############################################################################
#pca for temporal data weird
##############################################################################

cc<- as.data.frame(rownames(A1))
ph<- pheno[,c(3,27)]
dim(pca_result$x)
colnames(cc)<- 'GID'
ph<- ph %>% unique()

A2<- merge(cc, ph, by= 'GID')

s_df <- A2 %>% distinct(GID, .keep_all = TRUE)
pca_data$ped <-  s_df$Pedigree
# Display the resulting data frame with only unique values in the 'GID' column
print(unique_values_df)


# Specify the columns for PCA 
columns_for_pca <- spc21_2[, 5:8]

# Perform PCA
pca_result <- prcomp(columns_for_pca, scale. = TRUE)

# Create a data frame for the PCA results
pca_data <- as.data.frame(pca_result$x)

# Add 'JD' and 'Env' from the original data frame to the PCA data frame
pca_data$ped <- spc21_2$JD
pca_data$Env <- spc21_2$Env

# Create a PCA plot using the first two principal components
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = ped)) +
  geom_point() +
  labs(title = "PCA Plot 21_2") +
  theme_minimal()

# Print the PCA plot
print(pca_plot)
dim(pca_data)

pca_data<- pca_data %>% filter(!ped == c('Prop') )

####################################################################################################
#WEIRD PLOTTING ATTEMPTS OF THE NDVI-S INSIDE eNV AS LINES. 
####################################################################################################
ggplot(spc22, aes(x = NDVI, fill= JD)) +
  #facet_grid(~JD)+ 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_histogram(alpha = 0.5) +
  labs(fill = "JD")#+ 
  #theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #ggtitle("Prediction accuracy of NDVI phenomic predictions")



ggplot(spctr, aes(x= JD ,y = NDVI,  group= JD ,fill= JD)) +
  facet_grid(~Env)+ 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_boxplot() +
  labs(fill = "JD")#+ 
#theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#ggtitle("Predictio

ggplot(spctr, aes(x = NDVI,  group= JD ,fill= Env)) +
  facet_grid(~Env+ JD )+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_histogram() +
  labs(fill = "Env") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggtitle("Predictio

# Create a gradient color palette from light blue to dark blue
my_color_palette <- colorRampPalette(c("grey" , "black"))(7)

ggplot(spc22, aes(x = NDVI, group = JD, fill=  JD)) +
  facet_grid(~JD)+
  geom_density(alpha= 0.5) +
  #labs(fill = "JD") +
  scale_color_manual(values = my_color_palette)

#############################################################################################################################################
#Visualize top and bottom yield with NDVI
#############################################################################################################################################
NDlong<- MCGm %>%  pivot_longer(
  cols = starts_with("NDVI"),    # Columns to pivot
  names_to = "Timepoint",        # New column for time points
  values_to = "NDVI_value"       # New column for NDVI values
)

NDMlong<- MCGmM %>%  pivot_longer(
  cols = starts_with("NDVI"),    # Columns to pivot
  names_to = "Timepoint",        # New column for time points
  values_to = "NDVI_value"       # New column for NDVI values
)



NDlong$plot<- as.factor(NDlong$plot)
NDMlong$plot<- as.factor(NDMlong$plot)
color_palette <- scales::hue_pal()(nrow(NDVIMCG))


########################################################
# Sort the data by yield_kgha in descending order
sorted_data <- NDlong[order(-NDlong$yield_kgha), ]
# Select the top 10 highest and lowest rows
top10_highest <- head(sorted_data, 100)
top10_lowest <- tail(sorted_data, 100)
# Assign different colors to the subsets
top10_highest$color <- "High Yield"
top10_lowest$color <- "Low Yield"
# Combine the selected rows
selected_rows <- rbind(top10_highest, top10_lowest)
########################################

# Sort the data by yield_kgha in descending order
sorted_data <- NDMlong[order(-NDMlong$yield_kgha), ]
# Select the top 10 highest and lowest rows
top10_highest <- head(sorted_data, 100)
top10_lowest <- tail(sorted_data, 100)
# Assign different colors to the subsets
top10_highest$color <- "High Yield"
top10_lowest$color <- "Low Yield"

# Combine the selected rows
selected_rowsM <- rbind(top10_highest, top10_lowest)

# Create the plot with geom_smooth() lines
p1<- ggplot(selected_rows, aes(x = Timepoint, y = NDVI_value, group = plot, color = color)) +
  geom_line() +
  geom_smooth(size=0.8) +
  scale_color_manual(values = c("High Yield" = "blue", "Low Yield" = "red")) +
  theme_minimal()  +
  theme(panel.background = element_rect(fill = "gray")) +
  labs(title = "Non-masked plot NDVI High and Low Yielding Plots")



p2<- ggplot(selected_rowsM, aes(x = Timepoint, y = NDVI_value, group = plot, color = color)) +
  geom_line() +
  geom_smooth(size=0.8) +
  scale_color_manual(values = c("High Yield" = "blue", "Low Yield" = "red")) +
  theme_minimal()  +
  theme(panel.background = element_rect(fill = "gray")) +
  labs(title = "Masked plot NDVI High and Low Yielding Plots")

library(gridExtra)

grid.arrange(p1, p2, ncol=2)

ggplot(NDlong, aes(x=Timepoint, y=NDVI_value, group=plot)) +
  geom_smooth(aes(color= plot))+
  scale_color_manual(values = color_palette) +
  theme(legend.position = "none")


######################################################################################################################

