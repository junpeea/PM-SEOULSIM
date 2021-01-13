############################################################################################################
# RESEARCH THESIS: 
# Simulation studies between exposure to airborne particulate matter and low birth weight in Seoul, Korea, 2010.
# Type: Function code


# 20180101
# Yoon Bae Jun
############################################################################################################

############################################################################################################
# grfPM10 = function(predict.grid,pred_Z,variogpars,seed)
############################################################################################################
grfPM10 = function(predict.grid,GISdata,DongCenter,variogpars,seed){
  
  grid_Z1 <- as.numeric(gsub(",","",GISdata[,"Road_L_0100"]))
  grid_Z2 <- as.numeric(gsub(",","",GISdata[,"LS700_0500"]))
  grid_Z3 <- as.numeric(gsub(",","",GISdata[,"B_bnu_4_1000"]))
  grid_Z4 <- as.numeric(gsub(",","",GISdata[,"D_Bus"]))
  grid_Z5 <- as.numeric(gsub(",","",GISdata[,"B_bem_4_1000"]))
  grid_Z1 <- (grid_Z1 - mean(grid_Z1))/sd(grid_Z1)
  grid_Z2 <- (grid_Z2 - mean(grid_Z2))/sd(grid_Z2)
  grid_Z3 <- (grid_Z3 - mean(grid_Z3))/sd(grid_Z3)
  grid_Z4 <- (grid_Z4 - mean(grid_Z4))/sd(grid_Z4)
  grid_Z5 <- (grid_Z5 - mean(grid_Z5))/sd(grid_Z5)
  grid_Z  <- as.matrix(cbind(GISdata[,c("TM_X","TM_Y")],1,grid_Z1,grid_Z2,grid_Z3,grid_Z4,grid_Z5))
  
  dgrid <- DongCenter[,c("Road_L_0100","LS700_0500","B_bnu_4_1000","D_Bus","B_bem_4_1000")]
  refgrd <- (dgrid - matrix(rep(colMeans(dgrid),nrow(dgrid)),nrow=nrow(dgrid),byrow=T))/
    matrix(rep(apply(dgrid,2,sd),nrow(dgrid)),nrow=nrow(dgrid),byrow=T)
  dgrid_Z <- as.matrix(cbind(DongCenter[,c("TM_X","TM_Y")],1,refgrd))
  
  pred_Z = data.frame(rbind(grid_Z,as.matrix(cbind(rawTMSdata[,c("TM_X","TM_Y")],1,rawTMSdata[,7:11])),dgrid_Z))
  
  ind1 = paste0(pred_Z$TM_X,pred_Z$TM_Y)
  ind2 = paste0(predict.grid$TM.X,predict.grid$TM.Y)
  set.seed(seed)
  Mean = as.matrix(pred_Z[match(ind2,ind1),3:8])%*%matrix(variogpars[1:6],ncol=1)
  # Mean[which(is.na(Mean)==1)]
  grfoutput <- grf(nrow(predict.grid), grid=predict.grid[,1:2], mean=Mean, nugget=variogpars[7], cov.pars=variogpars[8:9], cov.m="exponential",messages=FALSE)
  return(grfoutput)
}

###########################################################################################################
# SimulateTMS = function(rawTMSdata,getPM10)
# PM10 Estimation Methods
#   SimpleAverage    = function(TMSdata,MyBirthdat)
#   NearestNeighbor  = function(TMSdata,MyBirthdat,GISdata)
#   IDWA             = function(TMSdata,MyBirthdat,GISdata)
#   LUR              = function(TMSdata,MyBirthdat,GISdata)
#   OrdinaryKriging  = function(TMSdata,MyBirthdat,GISdata)
#   UniversalKriging = function(TMSdata,MyBirthdat,GISdata)
# getEstPM10 = function(rawTMSdata,getPM10,MyBirthdat,GISdata) : True PM10 by the 6 prposed methods
############################################################################################################
SimulateTMS = function(rawTMSdata,getPM10){
  selectPMoverTMS = function(data.out,grf.out){
    ind01 = paste0(data.out[,"TM_X"],data.out[,"TM_Y"])
    ind02 = paste0(grf.out$coords[,1],grf.out$coords[,2])
    return(grf.out$data[match(ind01,ind02)])
  }
  dat = data.frame(cbind(selectPMoverTMS(rawTMSdata,getPM10[[1]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[2]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[3]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[4]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[5]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[6]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[7]]),
                         selectPMoverTMS(rawTMSdata,getPM10[[8]])))
  colnames(dat) <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                     "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  TMSdata = data.frame(cbind(rawTMSdata,dat))
  return(TMSdata)
}
# TMSdata = SimulateTMS(rawTMSdata,getPM10)

SimpleAverage = function(TMSdata,MyBirthdat,Messages=TRUE){
  n = nrow(MyBirthdat)
  colname <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
               "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  estPM10 <- matrix(0,n,8)
  for(i in 1:n){
    SGG <- MyBirthdat$SGG[i]
    ind1 = TMSdata$TMSflag==1
    ind2 = floor(TMSdata$EMD_cd/100)==SGG
    # ind2 = TMSdata$Sigungu_cd==SGG
    estPM10[i,] <- as.numeric(TMSdata[which(ind1 & ind2),colname])
    # estPM10[i,] <- as.numeric(colMeans(TMSdata[which(ind2),colname]))
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname
  return(estPM10)
}
# SimpleAverage(TMSdata,MyBirthdat,Messages=TRUE)

NearestNeighbor = function(TMSdata,MyBirthdat,GISdata,Messages=TRUE){
  n = nrow(MyBirthdat)
  colname <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
               "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  estPM10 <- matrix(0,n,8)
  for(i in 1:n){
    # ID <- MyBirthdat$GISID[i]
    ID <- match(myBirth2$GISID[i],GISdata$ID)
    loc <- c(GISdata$TM_X[ID],GISdata$TM_Y[ID])
    # ind1 = TMSdata$TMSflag==1
    # dist0 <- as.matrix(dist(rbind(TMSdata[ind1,c("TM_X","TM_Y")],loc)))
    dist0 <- as.matrix(dist(rbind(TMSdata[,c("TM_X","TM_Y")],loc)))
    distvec = dist0[ncol(dist0),1:(nrow(dist0)-1)]
    ind2 <- as.numeric(names(distvec[which(distvec==min(distvec))]))
    estPM10[i,] <- as.numeric(TMSdata[ind2,colname])
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname
  return(estPM10)
}

IDWA = function(TMSdata,MyBirthdat,GISdata,Messages=TRUE){
  n = nrow(MyBirthdat)
  colname <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
               "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  estPM10 <- matrix(0,n,8)
  for(i in 1:n){
    # ID <- MyBirthdat$GISID[i]
    ID <- match(myBirth2$GISID[i],GISdata$ID)
    loc <- c(GISdata$TM_X[ID],GISdata$TM_Y[ID])
    ind1 = TMSdata$TMSflag==1
    dist0 <- as.matrix(dist(rbind(TMSdata[ind1,c("TM_X","TM_Y")],loc)))
    distvec = dist0[ncol(dist0),1:(nrow(dist0)-1)]
    weight  = (1/distvec^2)/sum(1/distvec^2)
    estPM10[i,] <- as.numeric(t(weight)%*%as.matrix(TMSdata[ind1,colname]))
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname
  return(estPM10)
}

LUR = function(TMSdata,MyBirthdat,GISdata,Messages=TRUE){
  n = nrow(MyBirthdat)
  colname1 <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  colname2 <- c("model1_Par","model2_Par","model3_Par","model4_Par",
                "model5_Par","model6_Par","model7_Par","model8_Par")
  estPM10 <- matrix(0,n,8)
  fit1.LUR = glm(model1_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit2.LUR = glm(model2_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit3.LUR = glm(model3_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit4.LUR = glm(model4_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit5.LUR = glm(model5_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit6.LUR = glm(model6_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit7.LUR = glm(model7_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  fit8.LUR = glm(model8_True_PM10~geo1+geo2+geo3+geo4+geo5,family="gaussian",data=TMSdata)
  coef.LUR = cbind(coefficients(fit1.LUR),coefficients(fit2.LUR),
                   coefficients(fit3.LUR),coefficients(fit4.LUR),
                   coefficients(fit5.LUR),coefficients(fit6.LUR),
                   coefficients(fit7.LUR),coefficients(fit8.LUR))
  geoScal = function(geovec){
    output = list()
    output$mu = mean(geovec)
    output$sd = sd(geovec)
    output$scaled = (geovec-mean(geovec))/sd(geovec)
    return(output)
  }
  geo1 = geoScal(as.numeric(gsub(",","",GISdata$Road_L_010)))$scaled
  geo2 = geoScal(as.numeric(gsub(",","",GISdata$LS700_0500)))$scaled
  geo3 = geoScal(as.numeric(gsub(",","",GISdata$B_bnu_4_1000)))$scaled
  geo4 = geoScal(as.numeric(gsub(",","",GISdata$D_Bus)))$scaled
  geo5 = geoScal(as.numeric(gsub(",","",GISdata$B_bem_4_1000)))$scaled
  for(i in 1:n){
    # ID <- MyBirthdat$GISID[i]
    ID <- match(myBirth2$GISID[i],GISdata$ID)
    estPM10[i,] <- as.numeric(c(1,geo1[ID],geo2[ID],geo3[ID],geo4[ID],geo5[ID])%*%coef.LUR)
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname1
  estPara <- data.frame(cbind(c(coefficients(fit1.LUR),NA,NA,NA),c(coefficients(fit2.LUR),NA,NA,NA),
                              c(coefficients(fit3.LUR),NA,NA,NA),c(coefficients(fit4.LUR),NA,NA,NA),
                              c(coefficients(fit5.LUR),NA,NA,NA),c(coefficients(fit6.LUR),NA,NA,NA),
                              c(coefficients(fit7.LUR),NA,NA,NA),c(coefficients(fit8.LUR),NA,NA,NA)))
  colnames(estPara) <- colname2
  return(list(estPM10=estPM10,estPara=estPara))
}
# temp=LUR(TMSdata,MyBirthdat,GISdata)
# summary(temp$estPM10)
# t(temp$estPara)

OrdinaryKriging = function(TMSdata,MyBirthdat,GISdata,ErrFix=TRUE,Messages=TRUE){
  minpsill = 1
  minrange = 1000
  n = nrow(MyBirthdat)
  m = nrow(GISdata)
  colname1 <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  colname2 <- c("model1_Par","model2_Par","model3_Par","model4_Par",
                "model5_Par","model6_Par","model7_Par","model8_Par")
  estPM10  <- matrix(0,n,8)
  KrigPM10 <- matrix(0,m,8)
  estPara  <- matrix(0,9,8)
  ind1 = TMSdata$TMSflag==1
  for(model in 1:8){
    geod = as.geodata(TMSdata[ind1,],coords.col = c(5:6), data.col = (11+model), covar.col = c("geo1","geo2","geo3","geo4","geo5"))
    fit <- likfit(geod,trend=trend.spatial(trend ="cte", geodata = geod),cov.model="exponential",ini.cov.pars=c(1,1000),messages = FALSE)
    estPara[,model] <- as.numeric(c(fit$beta,NA,NA,NA,NA,NA,fit$nugget,fit$cov.pars))
    if(ErrFix){
      if(as.numeric(estPara[8,model]) < minpsill){estPara[8,model] <- minpsill} 
      if(estPara[9,model] < minrange){estPara[9,model] <- minrange}
    }
    krige = krige.control(type.krige = "ok",trend.d = "cte",cov.model = "exponential",
                          nugget=estPara[7,model],cov.pars=estPara[8:9,model])
    krige.OK = krige.conv(geodata=geod,locations=data.frame(cbind(GISdata$TM_X,GISdata$TM_Y)),
                          krige=krige)
    
    KrigPM10[,model] = krige.OK$predict 
  }
  
  for(i in 1:n){
    # ID <- MyBirthdat$GISID[i]
    ID <- match(myBirth2$GISID[i],GISdata$ID)
    estPM10[i,] <- KrigPM10[ID,]
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname1
  estPara <- data.frame(estPara)
  colnames(estPara) <- colname2
  return(list(estPM10=estPM10,estPara=estPara))
}
# temp=OrdinaryKriging(TMSdata,MyBirthdat,GISdata)
# summary(temp$estPM10)
# t(temp$estPara)

UniversalKriging = function(TMSdata,MyBirthdat,GISdata,ErrFix=TRUE,Messages=TRUE){
  minpsill = 1
  minrange = 1000
  n = nrow(MyBirthdat)
  m = nrow(GISdata)
  colname1 <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  colname2 <- c("model1_Par","model2_Par","model3_Par","model4_Par",
                "model5_Par","model6_Par","model7_Par","model8_Par")
  estPM10  <- matrix(0,n,8)
  KrigPM10 <- matrix(0,m,8)
  estPara  <- matrix(0,9,8)
  geoScal = function(geovec){
    output = list()
    output$mu = mean(geovec)
    output$sd = sd(geovec)
    output$scaled = (geovec-mean(geovec))/sd(geovec)
    return(output)
  }
  lgeo1 = geoScal(as.numeric(gsub(",","",GISdata$Road_L_010)))$scaled
  lgeo2 = geoScal(as.numeric(gsub(",","",GISdata$LS700_0500)))$scaled
  lgeo3 = geoScal(as.numeric(gsub(",","",GISdata$B_bnu_4_1000)))$scaled
  lgeo4 = geoScal(as.numeric(gsub(",","",GISdata$D_Bus)))$scaled
  lgeo5 = geoScal(as.numeric(gsub(",","",GISdata$B_bem_4_1000)))$scaled
  # geod_COA = as.geodata(cbind(1,GISdata,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(384:385), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  geod_COA = as.geodata(cbind(1,GISdata,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(8:9), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  for(model in 1:8){
    geod     = as.geodata(TMSdata,coords.col = c(5:6), data.col = (13+model), covar.col = c("geo1","geo2","geo3","geo4","geo5"))
    fit <- likfit(geod,trend=trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),cov.model="exponential",ini.cov.pars=c(1,1000),messages=FALSE)
    estPara[,model] <- as.numeric(c(fit$beta,fit$nugget,fit$cov.pars))
    if(ErrFix){
      if(as.numeric(estPara[8,model]) < minpsill){estPara[8,model] <- minpsill} 
      if(estPara[9,model] < minrange){estPara[9,model] <- minrange}
    }
    krige = krige.control(type.krige = "ok",
                          trend.d = trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),
                          trend.l = trend.spatial(trend=~1+lgeo1+lgeo2+lgeo3+lgeo4+lgeo5,geodata=geod_COA),
                          cov.model = "exponential",nugget=estPara[7,model],cov.pars=estPara[8:9,model])
    krige.UK = krige.conv(geodata=geod,locations=data.frame(cbind(GISdata$TM_X,GISdata$TM_Y)),
                          krige=krige)
    
    KrigPM10[,model] = krige.UK$predict 
  }
  
  for(i in 1:n){
    # ID <- MyBirthdat$GISID[i]
    ID <- match(myBirth2$GISID[i],GISdata$ID)
    estPM10[i,] <- KrigPM10[ID,]
    # if(Messages){
    #   if(i%%100==0) cat(i,date(),'\n')
    # }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname1
  estPara <- data.frame(estPara)
  colnames(estPara) <- colname2
  return(list(estPM10=estPM10,estPara=estPara))
}
# temp=UniversalKriging(TMSdata,MyBirthdat,GISdata)
# summary(temp$estPM10)
# t(temp$estPara)

EstimatePM10 = function(rawTMSdata,getPM10,MyBirthdat,GISdata,DongCenter,GuCenter,GridCenter,messages=FALSE){
  output = list()
  
  TMSdata = SimulateTMS(rawTMSdata,getPM10)
  
  output[[1]] = SimpleAverage(TMSdata,MyBirthdat,Messages=messages)
  
  output[[2]] = NearestNeighbor(TMSdata,MyBirthdat,GISdata,Messages=messages)
  
  output[[3]] = IDWA(TMSdata,MyBirthdat,GISdata,Messages=messages)
  
  output[[4]] = LUR(TMSdata,MyBirthdat,GISdata,Messages=messages)
  
  # output[[5]] = OrdinaryKriging(TMSdata,MyBirthdat,GISdata,Messages=messages)
  output[[5]] = UK.GCSimple(TMSdata,MyBirthdat,GuCenter,Messages=messages)
  
  output[[6]] = UniversalKriging(TMSdata,MyBirthdat,GISdata,Messages=messages)
  
  output[[7]] = UK.DCAverage(TMSdata,MyBirthdat,DongCenter,Messages=messages)
  
  output[[8]] = UK.COAverage(TMSdata,MyBirthdat,GISdata,Messages=messages)
  
  output[[9]] = UK.DCAverage(TMSdata,MyBirthdat,GridCenter,Messages=messages)
  
  return(output)
}

TrueSimPM10 = function(predict.grid,getPM10,MyBirthdat,Messages=TRUE){
  n = nrow(MyBirthdat)
  TM = matrix(0,n,8)
  # ind1 <- paste(predict.grid$TM.X,predict.grid$TM.Y)
  ind1 <- paste(predict.grid$TM_X,predict.grid$TM_Y)
  for(i in 1:n){
    # ID <- MyBirthdat$GISID[i]
    ID <- match(myBirth2$GISID[i],GISdata$ID)
    ind2 <- which(is.na(match(ind1,paste(GISdata$TM_X[ID],GISdata$TM_Y[ID])))==0)
    for(j in 1:8){
      TM[i,j] <- getPM10[[j]]$data[ind2]
    }
    if(Messages){
      if(i%%100==0) cat("TrueSimPM10     ",i,date(),'\n')
    }
  }
  return(TM)
}

# IFAddressBlinded = function(getEstPM10,SGGID){
#   
#   outEstPM10 <- getEstPM10
#   
#   outEstPM10[[1]] = getEstPM10[[1]]
#   
#   outEstPM10[[2]] <- getEstPM10[[2]] # No interpretation of this result!!!!!!!
#   
#   outEstPM10[[3]] <- getEstPM10[[3]] # No interpretation of this result!!!!!!!
#   
#   ExpectedSGG = function(get,SGGID){
#     out <- get
#     groupcount = table(SGGID)
#     groupnames = names(groupcount)
#     for(grp in 1:length(groupnames)){
#       ExpectedSGG = colMeans(get$estPM10[which(SGGID==groupnames[grp]),])
#       out$estPM10[which(SGGID==groupnames[grp]),] <-
#         matrix(rep(ExpectedSGG,groupcount[grp]),nrow=groupcount[grp],byrow=TRUE)
#     }
#     return(out)
#   }
#   
#   outEstPM10[[4]] = ExpectedSGG(getEstPM10[[4]],SGGID)
#   
#   outEstPM10[[5]] = ExpectedSGG(getEstPM10[[5]],SGGID)
#   
#   outEstPM10[[6]] = ExpectedSGG(getEstPM10[[6]],SGGID)
#   
#   return(outEstPM10)
# }

UK.GCSimple = function(TMSdata,MyBirthdat,GuCenter,ErrFix=TRUE,Messages=TRUE){
  minpsill = 1
  minrange = 1000
  n = nrow(MyBirthdat)
  m = nrow(GuCenter)
  colname1 <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  colname2 <- c("model1_Par","model2_Par","model3_Par","model4_Par",
                "model5_Par","model6_Par","model7_Par","model8_Par")
  estPM10  <- matrix(0,n,8)
  KrigPM10 <- matrix(0,m,8)
  estPara  <- matrix(0,9,8)
  geoScal = function(geovec){
    output = list()
    output$mu = mean(geovec)
    output$sd = sd(geovec)
    output$scaled = (geovec-mean(geovec))/sd(geovec)
    return(output)
  }
  lgeo1 = geoScal(as.numeric(gsub(",","",GuCenter$Road_L_010)))$scaled
  lgeo2 = geoScal(as.numeric(gsub(",","",GuCenter$LS700_0500)))$scaled
  lgeo3 = geoScal(as.numeric(gsub(",","",GuCenter$B_bnu_4_1000)))$scaled
  lgeo4 = geoScal(as.numeric(gsub(",","",GuCenter$D_Bus)))$scaled
  lgeo5 = geoScal(as.numeric(gsub(",","",GuCenter$B_bem_4_1000)))$scaled
  geod_GUC = as.geodata(cbind(1,GuCenter,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(320:321), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  for(model in 1:8){
    geod = as.geodata(TMSdata,coords.col = c(5:6), data.col = (13+model), covar.col = c("geo1","geo2","geo3","geo4","geo5"))
    fit <- likfit(geod,trend=trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),cov.model="exponential",ini.cov.pars=c(1,1000),messages=FALSE)
    estPara[,model] <- as.numeric(c(fit$beta,fit$nugget,fit$cov.pars))
    if(ErrFix){
      if(as.numeric(estPara[8,model]) < minpsill){estPara[8,model] <- minpsill} 
      if(estPara[9,model] < minrange){estPara[9,model] <- minrange}
    }
    krige = krige.control(type.krige = "ok",
                          trend.d = trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),
                          trend.l = trend.spatial(trend=~1+lgeo1+lgeo2+lgeo3+lgeo4+lgeo5,geodata=geod_GUC),
                          cov.model = "exponential",nugget=estPara[7,model],cov.pars=estPara[8:9,model])
    krige.UK = krige.conv(geodata=geod,locations=data.frame(cbind(GuCenter$TM_X,GuCenter$TM_Y)),
                          krige=krige)
    
    KrigPM10[,model] = krige.UK$predict 
  }
  for(i in 1:n){
    ID <- MyBirthdat$SGGID[i]
    estPM10[i,] <- KrigPM10[which(GuCenter$SGGID==ID),]
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname1
  estPara <- data.frame(estPara)
  colnames(estPara) <- colname2
  return(list(estPM10=estPM10,estPara=estPara))
}
# UK.GCSimple(TMSdata,MyBirthdat,GuCenter,Messages=TRUE)

UK.COAverage = function(TMSdata,MyBirthdat,GISdata,ErrFix=TRUE,Messages=TRUE){
  minpsill = 1
  minrange = 1000
  n = nrow(MyBirthdat)
  m = nrow(GISdata)
  colname1 <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  colname2 <- c("model1_Par","model2_Par","model3_Par","model4_Par",
                "model5_Par","model6_Par","model7_Par","model8_Par")
  estPM10  <- matrix(0,n,8)
  KrigPM10 <- matrix(0,m,8)
  estPara  <- matrix(0,9,8)
  geoScal = function(geovec){
    output = list()
    output$mu = mean(geovec)
    output$sd = sd(geovec)
    output$scaled = (geovec-mean(geovec))/sd(geovec)
    return(output)
  }
  lgeo1 = geoScal(as.numeric(gsub(",","",GISdata$Road_L_010)))$scaled
  lgeo2 = geoScal(as.numeric(gsub(",","",GISdata$LS700_0500)))$scaled
  lgeo3 = geoScal(as.numeric(gsub(",","",GISdata$B_bnu_4_1000)))$scaled
  lgeo4 = geoScal(as.numeric(gsub(",","",GISdata$D_Bus)))$scaled
  lgeo5 = geoScal(as.numeric(gsub(",","",GISdata$B_bem_4_1000)))$scaled
  # geod_COA = as.geodata(cbind(1,GISdata,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(384:385), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  geod_COA = as.geodata(cbind(1,GISdata,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(8:9), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  for(model in 1:8){
    geod = as.geodata(TMSdata,coords.col = c(5:6), data.col = (13+model), covar.col = c("geo1","geo2","geo3","geo4","geo5"))
    fit <- likfit(geod,trend=trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),cov.model="exponential",ini.cov.pars=c(1,1000),messages=FALSE)
    estPara[,model] <- as.numeric(c(fit$beta,fit$nugget,fit$cov.pars))
    if(ErrFix){
      if(as.numeric(estPara[8,model]) < minpsill){estPara[8,model] <- minpsill} 
      if(estPara[9,model] < minrange){estPara[9,model] <- minrange}
    }
    krige = krige.control(type.krige = "ok",
                          trend.d = trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),
                          trend.l = trend.spatial(trend=~1+lgeo1+lgeo2+lgeo3+lgeo4+lgeo5,geodata=geod_COA),
                          cov.model = "exponential",nugget=estPara[7,model],cov.pars=estPara[8:9,model])
    krige.UK = krige.conv(geodata=geod,locations=data.frame(cbind(GISdata$TM_X,GISdata$TM_Y)),
                          krige=krige)
    
    KrigPM10[,model] = krige.UK$predict 
  }
  for(i in 1:n){
    ID <- MyBirthdat$SGGID[i]
    # estPM10[i,] <- colMeans(KrigPM10[which(GISdata$SGG_cd_10==ID),])
    estPM10[i,] <- colMeans(KrigPM10[which(GISdata$SGGCD==ID),])
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname1
  estPara <- data.frame(estPara)
  colnames(estPara) <- colname2
  return(list(estPM10=estPM10,estPara=estPara))
}
# UK.COAverage(TMSdata,MyBirthdat,GISdata,ErrFix=TRUE,Messages=TRUE)

UK.DCAverage = function(TMSdata,MyBirthdat,DongCenter,ErrFix=TRUE,Messages=TRUE){
  minpsill = 1
  minrange = 1000
  n = nrow(MyBirthdat)
  m = nrow(DongCenter)
  colname1 <- c("model1_True_PM10","model2_True_PM10","model3_True_PM10","model4_True_PM10",
                "model5_True_PM10","model6_True_PM10","model7_True_PM10","model8_True_PM10")
  colname2 <- c("model1_Par","model2_Par","model3_Par","model4_Par",
                "model5_Par","model6_Par","model7_Par","model8_Par")
  
  crs <- '+proj=tmerc +lat_0=38 +lon_0=127.0028902777778 +k=1 +x_0=200000 +y_0=500000 +ellps=bessel +towgs84=-146.414,507.337,680.507,0,0,0,0 +units=m +no_defs'
  em.s <- SpatialPointsDataFrame(coords=DongCenter[,c('TM_X','TM_Y')],
                                 data=DongCenter,
                                 proj4string=CRS(crs))
  DongCenter$Sigungu_cd <- em.s %over% gu.poly.11
  
  estPM10  <- matrix(0,n,8)
  KrigPM10 <- matrix(0,m,8)
  estPara  <- matrix(0,9,8)
  geoScal = function(geovec){
    output = list()
    output$mu = mean(geovec)
    output$sd = sd(geovec)
    output$scaled = (geovec-mean(geovec))/sd(geovec)
    return(output)
  }
  lgeo1 = geoScal(as.numeric(gsub(",","",DongCenter$Road_L_010)))$scaled
  lgeo2 = geoScal(as.numeric(gsub(",","",DongCenter$LS700_0500)))$scaled
  lgeo3 = geoScal(as.numeric(gsub(",","",DongCenter$B_bnu_4_1000)))$scaled
  lgeo4 = geoScal(as.numeric(gsub(",","",DongCenter$D_Bus)))$scaled
  lgeo5 = geoScal(as.numeric(gsub(",","",DongCenter$B_bem_4_1000)))$scaled
  # geod_DC = as.geodata(cbind(1,DongCenter,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(11:12), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  geod_DC = as.geodata(cbind(1,DongCenter,lgeo1,lgeo2,lgeo3,lgeo4,lgeo5),coords.col = c(3:4), data.col = 1, covar.col = c("lgeo1","lgeo2","lgeo3","lgeo4","lgeo5"))
  for(model in 1:8){
    geod = as.geodata(TMSdata,coords.col = c(5:6), data.col = (13+model), covar.col = c("geo1","geo2","geo3","geo4","geo5"))
    fit <- likfit(geod,trend=trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),cov.model="exponential",ini.cov.pars=c(1,1000),messages=FALSE)
    estPara[,model] <- as.numeric(c(fit$beta,fit$nugget,fit$cov.pars))
    if(ErrFix){
      if(as.numeric(estPara[8,model]) < minpsill){estPara[8,model] <- minpsill} 
      if(estPara[9,model] < minrange){estPara[9,model] <- minrange}
    }
    krige = krige.control(type.krige = "ok",
                          trend.d = trend.spatial(trend=~1+geo1+geo2+geo3+geo4+geo5,geodata=geod),
                          trend.l = trend.spatial(trend=~1+lgeo1+lgeo2+lgeo3+lgeo4+lgeo5,geodata=geod_DC),
                          cov.model = "exponential",nugget=estPara[7,model],cov.pars=estPara[8:9,model])
    krige.UK = krige.conv(geodata=geod,locations=data.frame(cbind(DongCenter$TM_X,DongCenter$TM_Y)),
                          krige=krige)
    
    KrigPM10[,model] = krige.UK$predict 
  }
  for(i in 1:n){
    ID <- MyBirthdat$SGGID[i]
    estPM10[i,] <- colMeans(KrigPM10[which(DongCenter$Sigungu_cd==ID),])
    if(Messages){
      if(i%%100==0) cat(i,date(),'\n')
    }
  }
  
  estPM10 <- data.frame(estPM10)
  colnames(estPM10) <- colname1
  estPara <- data.frame(estPara)
  colnames(estPara) <- colname2
  return(list(estPM10=estPM10,estPara=estPara))
}
# UK.DCAverage(TMSdata,MyBirthdat,DongCenter,ErrFix=TRUE,Messages=TRUE)


