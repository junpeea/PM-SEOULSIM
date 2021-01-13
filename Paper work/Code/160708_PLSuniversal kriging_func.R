################################################################################
# 
# plsUK.func : universal kriging model fitting
# plsCV : pls cv
# cv.gr : cv group creation
# allsubset.reg : all subset regression (exhaustive search)
#
################################################################################

library(pls)
library(geoR)


##### Covariate data processing ################################################
covar.process <- function(covar.data = covar, exclude = TRUE, u.p =0.1, covar.list){

# log transformation & truncation
    dat2 <- covar.data
    I.var0 <- grepl("^D_", names(dat2))
    dat11 <- dat2[,I.var0]
    dat11[dat11 > 25000] <- 25000
    dat11[dat11 < 10] <- 10
    dat12 <- log(dat11)
    dat2[,I.var0] <- dat12          

    if(sum(grepl("EM",covar.list)) > 0){
        I.var1 <- grepl("^EM_", names(dat2))
        dat13 <- dat2[,I.var1]
        dat14 <- log(dat13)
        dat2[,I.var1] <- dat14          
    }
    var.list <- names(dat2)[grepl(paste(covar.list,collapse="|"), names(dat2))]       
    
    if(exclude){
# exclude variables with values less than 10% of being different from most common values 
        exclude1.var <- names(covar.data[,var.list])[apply(dat2[,var.list], 2, 
                          function(x) length(unique(x))/nrow(dat2) < u.p)]
# exclude variables max landuse variable < 10%
        LUR.list <- names(dat2)[grepl("^LS",names(dat2))]
        exclude2.var <- names(dat2[,LUR.list])[apply(dat2[,LUR.list], 2, 
                          function(x) max(x) < .1)]

        exclude.var <- c(exclude1.var, exclude2.var)
        var.list1 <- var.list[!(var.list %in% exclude.var)]
        covar.pro <- dat2[,!c(names(dat2) %in% exclude.var)]
    } else {
        covar.pro <- dat2
    }
    return(covar.pro)
}


##### pls CV ###################################################################
plsCV <- function(covar.data, exp.data, pol, n.comp = 5, covar.list,
              valid = "LOO", seg = 5){

    exp.data1 <- exp.data[match(covar.data$ID, exp.data$TMSID),]
    I.var <- grepl(covar.list, names(covar.data))
    covar.scale <- scale(as.matrix(covar.data[,I.var]))
    pls.data <- list(exp.data1[,pol], covar.scale[,!is.na(apply(covar.scale,2,sum))])
    names(pls.data) <- c('y','X')
    if(valid=="LOO"){
        plscv <- plsr(y ~ X, ncomp=n.comp, validation=valid, data=pls.data)
    } else {
        plscv <- plsr(y ~ X, ncomp=n.comp, validation=valid, segments=seg, data=pls.data)
    }
    return(plscv)
}


##### PLS regression ###########################################################
pls <- function(covar.data, covar.new = NULL, exp.data, pol, covar.list, 
             n.comp = 5, reest = NULL, cvGr = NULL){

    exp.data2 <- exp.data[match(covar.data$ID, exp.data$TMSID),]
    I.var <- grepl(covar.list, names(covar.data))
    covar.scale <- scale(as.matrix(covar.data[,I.var]))
    covar.m <- apply(covar.data[,I.var], 2, mean)
    covar.sd <- apply(covar.data[,I.var], 2, sd)
    covar.stat <- data.frame(Avg = covar.m, SD = covar.sd)

    if(is.null(covar.new)){
        covar.new <- covar.data
        new.data <- covar.scale
    } else {
        I.var <- grepl(covar.list, names(covar.new))
        X <- covar.new[,colnames(covar.scale)]
        m  <- matrix( rep(covar.m, each=nrow(X)), nrow=nrow(X) )
        s <- matrix( rep(covar.sd, each=nrow(X)), nrow=nrow(X) )
        new.data <- as.matrix((X - m)/s)
    }

    pls.out <- list()
    if(is.null(reest)){
        pls.data <- list(exp.data2[,pol], covar.scale)
        names(pls.data) <- c('y','X')
        pls <- plsr(y ~ X, ncomp=n.comp, data=pls.data)
        pls.pred <- predict(pls, type="scores", newdata=new.data)
    } else {
        pls.pred <- matrix(NA, nrow(exp.data2), n.comp)
        if(reest== "LOO"){
            cvind <- 1:nrow(exp.data2)        
        } else if(reest== "CV") {
            cvind <- cvGr
        }    
        for(i in unique(cvind)){
            i.train <- cvind!=i
            i.test <- cvind==i
            train <- list(exp.data2[i.train, pol], covar.scale[i.train,])
            test <- covar.scale[i.test,,drop=F]
            names(train) <- c('y','X')
            pls <- plsr(y ~ X, ncomp=n.comp, data=train)
            pls.pred[i.test,] <- predict(pls, type="scores", newdata=test)
        }
        colnames(pls.pred) <- paste(rep("Comp",n.comp),(1:n.comp),sep=".")
    }
    pred <- data.frame(covar.new[,!I.var], pls.pred)
    pls.out <- list(fit = pls, pred = pred)
    out <- list(covar.stat = covar.stat, pls.out = pls.out)
    return(out)
}


##### universal model (including pls) fitting ##################################
predCV <- function(cov.data, exp.data, pol, log.trans = FALSE, var.list, cvGr,
                  ncv = 5, ini.covp = c(10,10)){
# 1. data
    id <- cov.data$ID
    I.a <- cov.data$CODE=="A"
    I.b <- cov.data$CODE=="B"
    exp.data2 <- exp.data[match(id, exp.data$TMSID),]

    pred.cv1 = pred.cv1.parA = pred.cv2 = pred.cv2.parA = rep(NA, length(id))
    var.name1 <- paste(var.list, collapse="+")
    var.name2 <- paste(var.list, collapse="|")

# 3. LUR
    lur.par <- list()
    dat <- data.frame(cov.data, exp.data2[,pol,drop=F])
    par0 <- lm(as.formula(paste(pol,"~", var.name1)), data=dat)$coef
    for(i in 1:max(cvGr)){
        train <- dat[cvGr!=i,]
        test <- dat[cvGr==i,]
        train.lm <- lm(as.formula(paste(pol,"~", var.name1)), data=train)
        par.cv0 <- train.lm$coef

        train.A <- dat[cvGr!=i & I.a,]
        train.A.lm <- lm(as.formula(paste(pol,"~", var.name1)), data=train.A)
        par.cv0.A <- train.A.lm$coef
# pred          
        X <- cbind(1, as.matrix(test[,names(train)[grepl(var.name2,names(train))]]))
        pred.cv1[cvGr==i] <- X %*% par.cv0
        pred.cv1.parA[cvGr==i] <- X %*% par.cv0.A
        if(i==1) par.cv1 <- par.cv0 else par.cv1 <- rbind(par.cv1, par.cv0)
        if(i==1) par.cv1.A <- par.cv0.A else par.cv1.A <- rbind(par.cv1.A, par.cv0.A)
    }
    if(log.trans){
        dat[,pol] <- exp(dat[,pol])
        pred.cv1 <- exp(pred.cv1)
        pred.cv1.parA <- exp(pred.cv1.parA)
    }
    mse.cv1 <- mean( (dat[,pol] - pred.cv1)^2 ) 
    mse.cv1.parA  <- mean( (dat[,pol] - pred.cv1.parA)^2 ) 
    r2.cv1 <- 1 - mse.cv1/var(dat[,pol])
    r2.cv1.parA  <- 1 - mse.cv1.parA/var(dat[,pol])
    lur.cvstat <- c(mse.cv1, mse.cv1.parA, r2.cv1, r2.cv1.parA)
    par1 <- rbind(par0, par.cv1, par.cv1.A)
    rownames(par1) <- c("all", paste(rep(c("CV.parAll","CV.parA"),each=max(cvGr)), 
                        rep(paste(rep("cv",max(cvGr)),1:max(cvGr), sep=""), 2), sep=":"))
    colnames(par1) <- c("b0", var.list)
    lur.par <- par1
    lur.pred <- data.frame(exp.data2$TMSID, dat[,pol], pred.cv1)
    lur.pred.parA   <- data.frame(exp.data2$TMSID, dat[,pol], pred.cv1.parA)
    
# 4. UK
    uk.par <- list()
    uk.lik <- list()
    uk.geodata <- list()
    dat <- data.frame(cov.data, exp.data2[,pol,drop=F])
    dat$TM_X <- dat$TM_X/1000
    dat$TM_Y <- dat$TM_Y/1000
    coords.c <- grep("TM", names(dat))
    data.c <- grep(paste("^",pol,"$",sep=""), names(dat))
    covar.c <- grep(var.name2, names(dat))
    dat.geo <- as.geodata(dat, coords.col=coords.c, data.col=data.c, covar.col=covar.c)
    lik0 <- likfit(dat.geo, trend=trend.spatial(as.formula(paste("~",var.name1)), dat.geo), 
                   ini.cov.pars = ini.covp, cov.m="exp")
    par0 <- unlist(lik0[c("beta","phi","sigmasq","tausq")])
    ini.covp2 <- c(pmax(ini.covp[1], lik0$cov.pars[1]), pmax(ini.covp[2], lik0$cov.pars[2]))
    for(i in 1:max(cvGr)){
# par.est
        train <- dat[cvGr!=i,]
        test <- dat[cvGr==i,]
        train.geo <- as.geodata(train, coords.col=coords.c, data.col=data.c, covar.col=covar.c)
        test.geo <- as.geodata(test, coords.col=coords.c, data.col=data.c, covar.col=covar.c)
        lik <- likfit(train.geo, trend=trend.spatial(as.formula(paste("~",var.name1)), train.geo), 
                      ini.cov.pars= ini.covp2, cov.m="exp")

        train.A <- dat[cvGr!=i & I.a,]
        test.A <- dat[cvGr==i & I.a,]
        train.A.geo <- as.geodata(train.A, coords.col=coords.c, data.col=data.c, covar.col=covar.c)
        lik.A <- likfit(train.A.geo, trend=trend.spatial(as.formula(paste("~",var.name1)), train.A.geo), 
                        ini.cov.pars= ini.covp2, cov.m="exp")
# pred          
        pred.cv2[cvGr==i] <- krige.conv( train.geo, loc=test[,c("TM_X","TM_Y")], 
                         krige=krige.control(obj.m=lik, 
                         trend.d=trend.spatial(as.formula(paste("~",var.name1)), train.geo), 
                         trend.l=trend.spatial(as.formula(paste("~",var.name1)), test.geo)) )$predict
        pred.cv2.parA[cvGr==i] <- krige.conv( train.geo, loc=test[,c("TM_X","TM_Y")], 
                         krige=krige.control(obj.m=lik.A, 
                         trend.d=trend.spatial(as.formula(paste("~",var.name1)), train.geo), 
                         trend.l=trend.spatial(as.formula(paste("~",var.name1)), test.geo)) )$predict
        par.cv20 <- unlist(lik[c("beta","phi","sigmasq","tausq")])
        par.cv20.A <- unlist(lik.A[c("beta","phi","sigmasq","tausq")])
        if(i==1) par.cv2 <- par.cv20 else par.cv2 <- rbind(par.cv2, par.cv20)
        if(i==1) par.cv2.A <- par.cv20.A else par.cv2.A <- rbind(par.cv2.A, par.cv20.A)
    }
    if(log.trans){
        dat[,pol] <- exp(dat[,pol])
        pred.cv2 <- exp(pred.cv2)
        pred.cv2.parA <- exp(pred.cv2.parA)
    }
    mse.cv2 <- mean( (dat[,pol] - pred.cv2)^2 ) 
    mse.cv2.parA <- mean( (dat[,pol] - pred.cv2.parA)^2 ) 
    r2.cv2 <- 1 - mse.cv2/var(dat[,pol])
    r2.cv2.parA <- 1 - mse.cv2.parA/var(dat[,pol])

    uk.cvstat <- c(mse.cv2, mse.cv2.parA, r2.cv2, r2.cv2.parA)
    par2 <- rbind(par0, par.cv2, par.cv2.A)
    rownames(par2) <- c("all", paste(rep(c("CV.parAll","CV.parA"),each=max(cvGr)), 
                        rep(paste(rep("cv",max(cvGr)),1:max(cvGr), sep=""), 2), sep=":"))
    colnames(par2) <- c("b0", var.list,"phi","sigmasq","tausq")
    uk.par <- par2
    uk.lik <- lik0
    uk.geodata <- dat.geo
    uk.pred <- data.frame(exp.data2$TMSID, dat[,pol], pred.cv2)
    uk.pred.parA <- data.frame(exp.data2$TMSID, dat[,pol], pred.cv2.parA)

# 5. output
    names(lur.cvstat) = names(uk.cvstat) = c("MSE","MSE.parA","R2","R2.parA")
    a <- c("ID","obs","pred")
    names(lur.pred) = names(lur.pred.parA) = names(uk.pred) = names(uk.pred.parA) = a 

    list( cv.gr = cvGr, UK.geodata = uk.geodata, UK.lik = uk.lik,
      LUR.par = lur.par, LUR.stat = lur.cvstat, 
      LUR.pred = lur.pred, LUR.pred.parA = lur.pred.parA, 
      UK.par = uk.par, UK.stat = uk.cvstat, 
      UK.pred = uk.pred, UK.pred.parA = uk.pred.parA)
}


##### prediction at unknown points #############################################
pred <- function(covar.subject, pls = TRUE, cov.data, model.out,
                 log.trans = FALSE, message = TRUE){
    par.est  <- model.out$UK.lik
    model.data  <- model.out$UK.geodata
    I.dup <- (duplicated(paste(covar.subject$TM_X, covar.subject$TM_Y)))

### 2. pls score prediction
    if(pls){             
        pls.data <- cov.data 
        var.name <- rownames(pls.data$pls.out$fit$coefficients)
        covar.mon.stat <- pls.data$covar.stat
        pls.fit  <- pls.data$pls.out$fit

        ind.var <- match(var.name, names(covar.subject))
        ind.var <- ind.var[!is.na(ind.var)]
        X <- covar.subject[,ind.var]
        mon.m  <- matrix( rep(covar.mon.stat$Avg, each=nrow(X)), nrow=nrow(X) )
        mon.sd <- matrix( rep(covar.mon.stat$SD, each=nrow(X)), nrow=nrow(X) )
        X.scale <- as.matrix((X - mon.m)/mon.sd)
        cov.pred <- predict(pls.fit, type="scores", newdata=X.scale)
        n.comp <- ncol(pls.fit$scores)
        rownames(cov.pred) <- covar.subject$ID
        colnames(cov.pred) <- paste("Comp",1:n.comp,sep="")
### 2.1. variable selection
    } else {
        var.name <- cov.data$varL
        covar.mon.stat <- cov.data$covar.stat

        X <- covar.subject[,var.name,drop=F]   
        mon.m  <- matrix( rep(covar.mon.stat$Avg, each=nrow(X)), nrow=nrow(X) )
        mon.sd <- matrix( rep(covar.mon.stat$SD, each=nrow(X)), nrow=nrow(X) )
        cov.pred <- as.matrix((X - mon.m)/mon.sd)
        n.comp <- length(var.name)
        rownames(cov.pred) <- covar.subject$ID
    }    

### 3. exposure prediction
    colnames(model.data$covariate) <- colnames(cov.pred)
    xy <- cbind(covar.subject$TM_X/1000, covar.subject$TM_Y/1000)
    xy1 <- xy[!I.dup,]
    cov.pred1 <- cov.pred[!I.dup,,drop=F]
    n <- floor(nrow(xy1)/1000)
    equa <- as.formula(paste("~",paste(colnames(cov.pred), collapse="+")))
    if(n==0){
        sub.data0 <- cbind(xy1, cov.pred1)
        sub.data <- as.geodata(sub.data0, coords.col=1:2, covar.col=3:(2+n.comp))
        colnames(sub.data$covariate) <- colnames(cov.pred)
        pred.sub <- krige.conv( model.data, loc=xy1,
                      krige=krige.control(obj.m=par.est, 
                        trend.d=trend.spatial(equa, model.data),
                        trend.l=trend.spatial(equa, sub.data)) )
        pred1 <- pred.sub$predict
    } else {
        for(i in 1:(n+1)){
            if(i < n+1){
                loop <- (1000*(i-1)+1):(1000*i)
            } else {
                loop <- (1000*n+1):nrow(xy1)
            }    
            sub.data0 <- cbind(xy1[loop,], cov.pred1[loop,])
            sub.data <- as.geodata(sub.data0, coords.col=1:2, covar.col=3:(2+n.comp))
            colnames(sub.data$covariate) <- colnames(cov.pred)
            pred.sub <- krige.conv( model.data, loc=xy1[loop,],
                          krige=krige.control(obj.m=par.est, 
                            trend.d=trend.spatial(equa, model.data),
                            trend.l=trend.spatial(equa, sub.data)) )
            pred0 <- pred.sub$predict
            if(i==1) pred1 <- pred0 else pred1 <- c(pred1,pred0)
            if(message) cat(n,i,date(),"\n", file="time.txt", append=T)
        }    
    }
    pred <- rep(NA, nrow(covar.subject))
    pred[!I.dup] <- pred1
    for(i in which(I.dup)){
        dup.xy <- covar.subject[i, c("TM_X","TM_Y")]
        I.dup1 <- (xy1[,1] == rep(dup.xy[[1]]/1000,nrow(xy1))) & 
                    (xy1[,2] == rep(dup.xy[[2]]/1000,nrow(xy1))) 
        pred[i] <- pred1[I.dup1] 
    }
    if(log.trans) pred <- exp(pred)
    names(pred) <- covar.subject$ID 
    cov.pred.sub <- data.frame(ID=covar.subject$ID, X=xy[,1], Y=xy[,2], cov.pred)
    return(list(pred.subject = pred, cov.pred.subject = cov.pred.sub))
}


##### forward variable selection #############################################
forward.s <- function(covar.data, exp.data, pol, covar.list, n.var = 5){
    exp.data2 <- exp.data[match(covar.data$ID, exp.data$TMSID),]
    V <- grepl(covar.list, names(covar.data))
    exp.covar <- data.frame(covar.data[,V], exp.data2[,pol,drop=F])
    exp.ind <- which(names(exp.covar)==pol)
    v1 <- names(which.max(abs(cor(exp.covar))[pol,-exp.ind]))     # 1.1. the first variable with the highest cor
    #anova(lm(exp.covar[,pol] ~ exp.covar[,v1]))                  # 1.2. make sure that f-test with the first variable is significant: otherwise stop
    v <- v1
    for(j in 1:(n.var-1)){
        vr <- names(exp.covar)[!(names(exp.covar) %in% c(v,pol))]
        for(i in vr){
            fn <- formula(paste(pol, "~", paste(c(v, i), collapse="+")))
            f11 <- anova(lm(fn, data=exp.covar))[(j+1),"F value"]
            if(i==vr[1]) f12 <- f11 else f12 <- c(f12,f11)
        }
        i1 <- abs(cor(exp.covar)[v,vr,drop=F]) < 0.7      # additonal step: do not select the variable highly correated with the first variable (=> 0.7)
        i2 <- apply(i1,2,sum) == length(v)
        v2 <- vr[i2][which.max(f12[i2])]                  # 2.1. the second variable with the highest partial f-test stat in the model with the first variable
        v <- c(v, v2)
    }
    return(v)   
}



#### cv group creation #########################################################
cv.gr <- function(LUR.data, n.cv){
  Ind <- 1:nrow(LUR.data)
  Ind.s <- Ind[sample.int(length(Ind),length(Ind))]
  Ind.max <- matrix(c(Ind.s,rep(NA,ceiling(length(Ind.s)/n.cv)*n.cv - length(Ind.s))), nrow=n.cv)
  Ind.cv <- matrix(FALSE, dim(LUR.data)[[1]],n.cv)
  for(i in 1:n.cv) {
    Ind.cv[,i] <- Ind %in% Ind.max[i,]
  }
  cv.LUR <- apply(Ind.cv,1,which)
}



#### exhaustive search #########################################################
allsubset.reg <- function( LUR.list, n.max, DATA, Y, cv){
  for(k in 1:n.max){
    lur.com <- combn(LUR.list, k)
    for(j in 1:ncol(lur.com)){
      lur.list2 <- lur.com[,j]
      dat1 <- DATA[,c(Y, lur.list2)]
      pred.cv <- rep(NA, nrow(DATA))
      for(i in 1:max(cv)){
        train <- dat1[cv!=i,]
        test <- dat1[cv==i,]
        train.lm <- lm(as.formula(paste(Y, "~", paste(lur.list2, collapse="+"))), data=train)
        X <- cbind(1,matrix(unlist(test[,lur.list2]), ncol=length(lur.list2)))
        pred.cv[cv==i] <- X %*% train.lm$coef
      }
      mse.cv <- mean( (DATA[,Y] - pred.cv)^2 ) 
      r2.cv <- 1 - mse.cv/var(DATA[,Y])
      cv.stat0 <- data.frame(k, paste(lur.list2, collapse=","), mse.cv, r2.cv)
      if(j==1){ 
        cv.stat1 <- cv.stat0
        pred1 <- pred.cv
      } else {
        cv.stat1 <- rbind(cv.stat1,cv.stat0)
        pred1 <- rbind(pred1,pred.cv)
      }  
    }  
    if(k==1){ 
      cv.stat <- cv.stat1
      pred <- pred1
    } else {
      cv.stat <- rbind(cv.stat,cv.stat1)
      pred <- rbind(pred, pred1)
    }   
  }
  names(cv.stat) <- c("var.n","var.list","MSE","R2")
  out <- list(predCV = pred, statCV = cv.stat)
  return(out)
}
