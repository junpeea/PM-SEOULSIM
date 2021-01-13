############################################################################################################
# RESEARCH THEME: 
# The impact of using a few fixed monitoring sites and partially known residential address
# on estimating health effects of PM10 on low birth weight in Seoul, Korea

# Source name : PM10_LBW_Main
# Author: Yoon Bae Jun/ Seoun Young Kim
# Date: 20201001
############################################################################################################

# Request: Environment settings
rm(list = ls())
require(geoR);require(grDevices);require(boot)
require(maptools);require(fields);require(maptools)
require(dplyr)
setwd("D:/KSY_myFirstPaper")
# source(file="PM10_LBW_func_180101.R")
# source(file="PM10_LBW_func_201001.R")
source(file="PM10_LBW_func_201101.R")

### Previously defined functions for prediction
source("160708_PLSuniversal kriging_func.R")

### load modeling outputs of PLS universal kriging based on monitoring data
load("ModelOutput.Rdata")

### Load the geographic variables
load('Song_and_Kim_2016_Data.RData')

# Load data-step
load(file="myLBW.Rdata")
load(file="myBirthdat.Rdata")
# load(file="GIS.Rdata")
load(file="rawTMS.Rdata")
# load(file="predict_grid.Rdata")
load(file="variogpars.Rdata") # (Reference: PM10_LBW_170830.R)
load(file="rawDC.Rdata")
load(file="Song_and_Kim_2016_Data.Rdata") # (Reference: Song and Kim 2016)
GISdata    = read.csv(file="SeoulResidCOA.csv")
GuCenter   = read.csv(file="gu_pos.csv")
# DongCenter2 = read.csv(file="SeoulGrid_Variable_150503.csv")
DongCenter = gv.em2
GridCenter = gv.em3
# temp = matrix(NA,length(which(EMD@data$SID_CD==11)),2)
# for (i in which(EMD@data$SID_CD==11)){
#   temp[i,] <- EMD@polygons[[i]]@Polygons[[1]]@coords %>% colMeans
# }
# summary(temp)
#################################### Korea shapefiles (Reference: KSY Project team)
# SIDO   <- readShapePoly("C:/Users/user/Desktop/Boundary_2010/SIDO_2010.shp")
SGG    <- readShapePoly("C:/Users/user/Desktop/Boundary_2010/BND_SIGUNGU_PG.shp")
EMD      <- readShapePoly("C:/Users/user/Desktop/Boundary_2010/EMD_2010.shp") #median(EMD$Shape_Area)/(1000^2) [1] 8.685693
# COA    <- readShapePoly("C:/Users/user/Desktop/Boundary_2010/COA_2010.shp")
# nrow(EMD@data[which(EMD$SID_CD==11),]) # [1] 424

# require(splines)
# formula = as.formula("lbw~pm10_oneyr+ns(ym,11)+as.factor(sex)+as.factor(m_edu)+as.factor(mage)+as.factor(parity)+as.factor(season)+ges_pr")
# formula = as.formula("lbw~pm10_oneyr")
# temp = summary(glm(formula,family=binomial,data=mydata))
# exp(IQR(mydata$pm10_oneyr)*temp$coefficients[2,1])
# exp(IQR(mydata$pm10_oneyr)*(temp$coefficients[2,1] - qnorm(0.975)*temp$coefficients[2,2]))
# exp(IQR(mydata$pm10_oneyr)*(temp$coefficients[2,1] + qnorm(0.975)*temp$coefficients[2,2]))


### True Health parameters
# LBW proportion (data-driven)
lbw.prop = length(which(mydata$lbw==1))/nrow(mydata) # [1] 0.01561102
# PM10 coefficient (Reference : Choe et al. 2018)
truebeta.coef = log(1.04)/IQR(mydata$pm10_oneyr) # 0.003275202
# intercept : Not available... -> arbitrary choosed with respect to available references...
truebeta.intercept = logit(lbw.prop) - truebeta.coef*50 # -4.307797
# Ture beta (~= log odds ratio ~= log relative risk)
truebeta = c(truebeta.intercept,truebeta.coef);truebeta # [1] -4.307797069  0.003275202

###########################################################################################################
# Stratified Random Sampling (SRS): 
############################################################################################################
myBirth2 = mydata %>% filter(b_yr==2010, mjob==0) %>% 
  select(obs,gu,sex,b_yr,b_mth,d_job,m_job,d_edu,m_edu,ges_pr,order,wgt,number,m_age,lbw,season,ges,mage,
         smoke,parity,brate,density,mjob,ym,income,bdate,fulldate,oneyr,trm2date,trm3date)
myBirth2$SGGID = myBirth2$gu + 11000
# myBirth2$EMDID = GISdata$ADM_DR_CD[match(myBirth2$GISID,as.numeric(GISdata$OBJECTID))]
# myBirth2$EMDID = GISdata$ADM_DR_CD[match(myBirth2$SGGID,as.numeric(GISdata$OBJECTID))]
n.birth = nrow(myBirth2)
seed <- 123

for(i in 1:n.birth){
  # sub.grid = GISdata[which(is.na(match(GISdata$SGG_cd_10,myBirth2$SGGID[i]))==0),]
  sub.grid = GISdata[which(is.na(match(GISdata$SGGCD,myBirth2$SGGID[i]))==0),]
  m <- nrow(sub.grid)
  p <- sub.grid$POP_1000/sum(sub.grid$POP_1000)
  set.seed(seed)
  ind <- sample.int(m,size=1,prob=p)
  seed <- seed + 10
  # myBirth2$GISID[i] <- sub.grid$OBJECTID[ind]
  myBirth2$GISID[i] <- sub.grid$ID[ind]
  if(i%%100==0) cat(i,date(),'\n')
}
table(myBirth2$GISID)

############################################################################################################
# PM10 Exposure - Health Effect Simulation
AddressBlind = FALSE; SGGID = myBirth2$SGGID; n = nrow(myBirth2)
# simnum = 1000; initialseed = 123
simnum = 1000; initialseed = 123
output.table2 = list(); output.table3 = list()
rawTMSdata$SGGCD = substr(rawTMSdata$EMD_cd,1,5)
rawTMSdata$ID    = rep(NA,nrow(rawTMSdata))
predict.grid2 = 
  rbind(
    GISdata[match(myBirth2$GISID,GISdata$ID),c("TM_X","TM_Y","SGGCD","ID")],
    rawTMSdata[,c("TM_X","TM_Y","SGGCD","ID")]
  )
colnames(predict.grid2) <- c("TM.X","TM.Y","SGGCD","ID")

############################################################################################################
start = Sys.time()
for(sim in 1:simnum){
  seed = initialseed + sim*100
  if(!file.exists(paste0("SimRdata/New_KSY_Simulation_seed(",seed,").Rdata"))){
    set.seed(seed)
    getPM10 = list()
    for(k in 1:8){
      getPM10[[k]] = grfPM10(predict.grid2,GISdata,DongCenter,variogpars[[k]],seed=seed)
      
    }
    getEstPM10  = EstimatePM10(rawTMSdata,getPM10,myBirth2,GISdata,DongCenter,GuCenter,GridCenter)
    # getTruePM10 = TrueSimPM10(predict.grid2,getPM10,myBirth2)
    # save(getPM10,getEstPM10,getTruePM10,file=paste0("SimRdata/New_KSY_Simulation_seed(",seed,").Rdata"))
    save(getPM10,getEstPM10,file=paste0("SimRdata/New_KSY_Simulation_seed(",seed,").Rdata"))
  }
}
record = Sys.time() - start
record

start = Sys.time()
for(sim in 1:simnum){
  seed = initialseed + sim*100
  if(file.exists(paste0("SimRdata/New_KSY_Simulation_seed(",seed,").Rdata"))){
    
    load(file=paste0("SimRdata/New_KSY_Simulation_seed(",seed,").Rdata"))
    
    getTruePM10 = data.frame(cbind(
      getPM10[[1]]$data,getPM10[[2]]$data,getPM10[[3]]$data,getPM10[[4]]$data,getPM10[[5]]$data,getPM10[[6]]$data,getPM10[[7]]$data,getPM10[[8]]$data))
    pmat = inv.logit(truebeta[1]+truebeta[2]*(getTruePM10 %>% as.matrix))
    model.lbw = matrix(0,n,8)
    for(i in 1:n){
      model.lbw[i,] <- c(rbinom(1,1,pmat[i,1]),rbinom(1,1,pmat[i,2]),rbinom(1,1,pmat[i,3]),rbinom(1,1,pmat[i,4]),
                         rbinom(1,1,pmat[i,5]),rbinom(1,1,pmat[i,6]),rbinom(1,1,pmat[i,7]),rbinom(1,1,pmat[i,8]))
    }
    
    TM  = getTruePM10[1:nrow(myBirth2),]
    # if(AddressBlind) getEstPM10 <- IFAddressBlinded(getEstPM10,SGGID)
    EM1 = getEstPM10[[1]]; dim(EM1)
    EM2 = getEstPM10[[2]]
    EM3 = getEstPM10[[3]]
    EM4 = getEstPM10[[4]]$estPM10
    EM5 = getEstPM10[[5]]$estPM10
    EM6 = getEstPM10[[6]]$estPM10
    EM7 = getEstPM10[[7]]$estPM10
    EM8 = getEstPM10[[8]]$estPM10
    EM9 = getEstPM10[[9]]$estPM10
    
    # output.Est  = matrix(NA,8*(8+1),6)
    output.Est  = matrix(NA,8*(9+1),6)
    
    for(model in 1:8){
      Model.Healthdata <- data.frame(cbind(model.lbw[,model],TM[,model],EM2[,model],EM3[,model],EM4[,model],EM6[,model],EM5[,model],EM1[,model],EM7[,model],EM8[,model],EM9[,model]))
      colnames(Model.Healthdata) <- c("lbw","TM","EM2","EM3","EM4","EM6","EM5","EM1","EM7","EM8","EM9")
      write.csv(Model.Healthdata,paste0(Sys.Date(),"_Healthdata_ES",model,".csv"))
    }
    
    for(model in 1:8){
      # Model.Healthdata <- data.frame(cbind(model.lbw[,model],TM[,model],EM2[,model],EM3[,model],EM4[,model],EM6[,model],EM5[,model],EM1[,model],EM7[,model],EM8[,model]))
      # colnames(Model.Healthdata) <- c("lbw","TM","EM2","EM3","EM4","EM6","EM5","EM1","EM7","EM8")
      Model.Healthdata <- data.frame(cbind(model.lbw[,model],TM[,model],EM2[,model],EM3[,model],EM4[,model],EM6[,model],EM5[,model],EM1[,model],EM7[,model],EM8[,model],EM9[,model]))
      colnames(Model.Healthdata) <- c("lbw","TM","EM2","EM3","EM4","EM6","EM5","EM1","EM7","EM8","EM9")
      fit = glm("lbw~1+TM",family="binomial",data=Model.Healthdata)
      # if(as.integer(summary(fit)$coefficients[2,4] < 0.10)){
      if(1){
        output.Est[(9+1)*model,1] = as.matrix(summary(fit)$coefficients)[2,1]
        output.Est[(9+1)*model,2] = as.matrix(summary(fit)$coefficients)[2,1]-qnorm(0.975)*as.matrix(summary(fit)$coefficients)[2,2]
        output.Est[(9+1)*model,3] = as.matrix(summary(fit)$coefficients)[2,1]+qnorm(0.975)*as.matrix(summary(fit)$coefficients)[2,2]
        output.Est[(9+1)*model,4] = (as.matrix(summary(fit)$coefficients)[2,1]-truebeta[2])^2 + as.matrix(summary(fit)$coefficients)[2,2]^2
        flag1 = as.integer(output.Est[(9+1)*model,2]<truebeta[2])
        flag2 = as.integer(truebeta[2]<output.Est[(9+1)*model,3])
        output.Est[(9+1)*model,5] = flag1*flag2
        output.Est[(9+1)*model,6] = as.integer(summary(fit)$coefficients[2,4] < 0.05)
        # if(output.Est[9*model,6]==1){
        # for(method in c(2,3,4,6,5,1,7,8)){
        for(method in c(2,3,4,6,5,1,7,8,9)){
          formula = as.formula(paste0("lbw~1+EM",method))
          fit = glm(formula,family="binomial",data=Model.Healthdata)
          # if(abs(as.matrix(summary(fit)$coefficients)[2,1])>1){
          #
          # }
          if(nrow(as.matrix(summary(fit)$coefficients))==2){
            output.Est[(9+1)*(model-1)+method,1] = as.matrix(summary(fit)$coefficients)[2,1]
            output.Est[(9+1)*(model-1)+method,2] = as.matrix(summary(fit)$coefficients)[2,1]-qnorm(0.975)*as.matrix(summary(fit)$coefficients)[2,2]
            output.Est[(9+1)*(model-1)+method,3] = as.matrix(summary(fit)$coefficients)[2,1]+qnorm(0.975)*as.matrix(summary(fit)$coefficients)[2,2]
            output.Est[(9+1)*(model-1)+method,4] = (as.matrix(summary(fit)$coefficients)[2,1]-truebeta[2])^2 + as.matrix(summary(fit)$coefficients)[2,2]^2
            flag1 = as.integer(output.Est[(9+1)*(model-1)+method,2]<truebeta[2])
            flag2 = as.integer(truebeta[2]<output.Est[(9+1)*(model-1)+method,3])
            output.Est[(9+1)*(model-1)+method,5] = flag1*flag2
            output.Est[(9+1)*(model-1)+method,6] = as.integer(summary(fit)$coefficients[2,4] < 0.05)
            if(abs(output.Est[(9+1)*(model-1)+method,1])>1){
              print(paste("sim:",sim,"model:",model,"method:",method,"beta:",summary(fit)$coefficients[2,1],"pvalue:",round(summary(fit)$coefficients[2,4],2)))
            }
            if(is.na(output.Est[(9+1)*(model-1)+method,1])==1){
              print(paste("sim:",sim,"model:",model,"method:",method,"beta:",summary(fit)$coefficients[2,1],"pvalue:",round(summary(fit)$coefficients[2,4],2)))
            }
            # if(round(summary(fit)$coefficients[2,4],2)>0.05){
            #   print(paste("sim:",sim,"model:",model,"method:",method,"beta:",summary(fit)$coefficients[2,1],"pvalue:",round(summary(fit)$coefficients[2,4],2)))
            # }
          }
        }
      }
    }

    # output.table2[[sim]] <- output.Est
    save(output.Est,paste0("SimRdata/New_outputEst_seed(",seed,").Rdata"))
    output.Pars = list()
    for(model in 1:8){
      output.Pars[[model]] = data.frame(rbind(t(getEstPM10[[4]]$estPara[,model]),
                                              t(getEstPM10[[5]]$estPara[,model]),
                                              t(getEstPM10[[6]]$estPara[,model]),
                                              t(getEstPM10[[7]]$estPara[,model]),
                                              t(getEstPM10[[8]]$estPara[,model]),
                                              t(getEstPM10[[9]]$estPara[,model])))
      rownames(output.Pars[[model]]) <- c("Method4","Method5","Method6","DA","CA","GridA")
      
    }
    # output.table3[[sim]] <- output.Pars
    save(output.Pars,paste0("SimRdata/New_outputPar_seed(",seed,").Rdata"))
    
  }
  
  cat('sim=',sim,' ',date(),'\n')
}

record = Sys.time() - start
record


for(sim in 1:simnum){
  seed = initialseed + sim*100
  load(Est,paste0("SimRdata/New_outputEst_seed(",seed,").Rdata"));  output.table2[[sim]] <- output.Est
  load(Est,paste0("SimRdata/New_outputPar_seed(",seed,").Rdata"));  output.table3[[sim]] <- output.Pars
}

############################################################ error correction procedures
output.corrected = output.table2
for(sim in 1:simnum){
  for(i in 1:28){
    if(is.na(abs(output.table2[[sim]][i,1]))==0){
      if(abs(output.table2[[sim]][i,1])>1) output.corrected[[sim]][i,] <- NA  
    }
  }
}
output.corrected[[sim]]

##################################################################### Finaloutput.table2

finaloutput.table2 = matrix(NA,8*9,6)
EstiMean = matrix(NA,8*9,simnum)
LowMean  = matrix(NA,8*9,simnum)
UppMean  = matrix(NA,8*9,simnum)
SSE      = matrix(NA,8*9,simnum)
Cpcount  = matrix(NA,8*9,simnum)
Spcount  = matrix(NA,8*9,simnum)

# abnorm.result=matrix(NA,simnum,2)
# norm.result  =matrix(NA,simnum,2)
# psillthresh = 1.0
for(sim in 1:simnum){
  # for(scene in 1:8){
  #   if(output.table3[[sim]][[scene]][2,8]<psillthresh){
  #     # abnorm.result[sim,1] <- output.corrected[[sim]][47,1]
  #     # abnorm.result[sim,2] <- output.table3[[sim]][[7]][2,8]
  #     output.corrected[[sim]][(7*scene-2),] <- rep(NA,5)
  #   }else{
  #     # norm.result[sim,1] <- output.corrected[[sim]][47,1]
  #     # norm.result[sim,2] <- output.table3[[sim]][[7]][2,8]
  #   }
  # }
  # if(output.corrected[[sim]][72,2]>0){
  if(1){
    
    EstiMean[,sim] <- output.corrected[[sim]][,1]
    LowMean[,sim]  <- output.corrected[[sim]][,2]
    UppMean[,sim]  <- output.corrected[[sim]][,3]
    SSE[,sim]      <- output.corrected[[sim]][,4]
    Cpcount[,sim]  <- output.corrected[[sim]][,5]
    Spcount[,sim]  <- output.corrected[[sim]][,6]
    
  }
  
}
# output.corrected[[2]][,2]
# summary(rowMeans(LowMean))
# finaloutput.table2[,1] = paste0(round(exp(rowMeans(EstiMean,na.rm=T)),2),"[",round(exp(rowMeans(LowMean,na.rm=T)),2),",",round(exp(rowMeans(UppMean,na.rm=T)),2),"]")
# finaloutput.table2[,1] = paste0("[",round(rowMeans(LowMean,na.rm=T),3),",",round(rowMeans(UppMean,na.rm=T),3),"]")
finaloutput.table2[,1] = round(( rowMeans(UppMean,na.rm=T) - rowMeans(LowMean,na.rm=T) ) / 2, 3)
finaloutput.table2[,2] = round(matrix(apply(EstiMean,1,function(w){return(sd(w,na.rm=T))}),ncol=1),3)
finaloutput.table2[,3] = round((rowMeans(EstiMean,na.rm=T) - truebeta[2])*100,2)
finaloutput.table2[,4] = round(rowMeans(SSE,na.rm=T)*100,3)
# for(scene in 1:8){
#   if(AddressBlind){
#     finaloutput.table2[(7*(scene-1)+1):(7*scene),5] <- round(c(
#       NA,NA,NA,
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+4,]))[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+5,]))[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+6,]))[1,2],NA
#     ),4)
#   }else{
#     finaloutput.table2[(7*(scene-1)+1):(7*scene),5] <- round(c(
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+1,]))[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+2,]))[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+3,]))[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+4,]))[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+5,]),use="complete.obs")[1,2],
#       cor(cbind(EstiMean[7*scene,],EstiMean[7*(scene-1)+6,]))[1,2],NA
#     ),4)
#   }
#   
# }
finaloutput.table2[,5] = round(rowMeans(Cpcount,na.rm=T),2)
finaloutput.table2[,6] = round(rowMeans(Spcount,na.rm=T),2)
finaloutput.table2
# write.csv(finaloutput.table2,"final_table2.csv")

finaloutput.table4 = matrix(NA,nrow=11*8,ncol=6)
for(i in 1:8){
  # finaloutput.table4[(11*(i-1)+1):(11*i),1:6] = finaloutput.table2[c(9,NA,2,3,4,6,NA,1,7,8,NA)+(9*(i-1)),]
  finaloutput.table4[(11*(i-1)+1):(11*i),1:6] = finaloutput.table2[c(9,NA,2,3,4,6,NA,1,5,7,8)+(9*(i-1)),]
}
finaloutput.table4
write.csv(finaloutput.table4,paste0(Sys.Date(),"_final_table4.csv"))
# write.csv(finaloutput.table4,"final_table5.csv")

Table4 = matrix(NA,8,10)
for(es in 1:8){
  temdat = data.frame(TM[1:46007,es],EM1[,es],EM2[,es],EM3[,es],EM4[,es],EM5[,es],EM6[,es],EM7[,es],EM8[,es],EM9[,es])
  Table4[es,] <- paste0(apply(temdat,2,mean) %>% round(2),"(",apply(temdat,2,sd) %>% round(2),")")
}
Table4
write.csv(Table4,paste0(Sys.Date(),"Table4.csv"))



