setwd("your_CD")
#install.packages('Matching')
#install.packages('survey')
#install.packages('MatchIt')
#install.packages('optmatch')

#library(Matching)
#library(survey)
#library(MatchIt)
#library(optmatch)

rm(list=ls())

source("~/R/sim_functions")

scenario.1 <- r.matrix(.3,.1,.1,.1,.1,.1,.1,.1,.1,.1)
##########
b0t <- 0.277777778
b0c <- 0
g <- 0.277777778
b1 <- 0.055555556
b2 <- 0.055555556
b3 <- 0.055555556
b4 <- 0.055555556

tau_list = c(0,2,3,4,5)
sigma.e=1

nT<-250
nC<-750
xC_bar<-0
xT_bar<-0
xC_sd<-1
xT_sd<-1

#sigma.y=sigma.g+sigma.e
nTreatClust=20
nControlClust=20
nChildPerClust=25

nClust <- nTreatClust + nControlClust
nObs <- nChildPerClust * nClust

S= 10000#number of replications

pb <- txtProgressBar(min = 0, max = S, style = 3)
#begin loop for tau_ICC
for (tau in tau_list){
  
  ICC<-matrix(NA, nrow=S,ncol=5) 
  colnames(ICC)<-c("sigma2y","sigma2alpha","icc","Fvalue","pvalue")
  
  Mtch<-matrix(NA, nrow=S, ncol=12)
  colnames(Mtch)<-c("ate","sd.ate","ols","sd.ols","psm.ate","sd.psm","psm.att","sd.att","ipw",
                    "sd.ipw","sub.ate","sd.sub")
  
  #begin simulation of data sets 
  
  for(s in 1:S) {  
    set.seed(s)
    #generate variables
    x1_t <- rnorm(nT, mean=xT_bar, sd=xT_sd)
    x2_t <- rnorm(nT, mean=xT_bar, sd=xT_sd)
    x3_t <- rnorm(nT, mean=xT_bar, sd=xT_sd)
    x4_t <- rnorm(nT, mean=xT_bar, sd=xT_sd)
    
    x1_c <- rnorm(nC, mean=xC_bar, sd=xC_sd)
    x2_c <- rnorm(nC, mean=xC_bar, sd=xC_sd)
    x3_c <- rnorm(nC, mean=xC_bar, sd=xC_sd)
    x4_c <- rnorm(nC, mean=xC_bar, sd=xC_sd)
    
    t <-c(rep(1, nT), rep(0, nC))  # make trt indicator #
    
    x1 <-(c(x1_t,x1_c))
    x2 <-(c(x2_t,x2_c))
    x3 <-(c(x3_t,x3_c))
    x4 <-(c(x4_t,x4_c))
    #assing desired correlations
    d <- matrix(c(t,x1,x2,x3,x4), nrow=nObs, ncol=5)
    colnames(d)<-c("tx","x1","x2","x3","x4")
    #sim.d[,1] <- ifelse(sim.d[,1] > mean(sim.d[,1]), 1, 0)
    #x[,1] <- cut.v(x[,1], c(0, 1-p$treat, 1)) - 1
    x.cor <- d %*% chol(scenario.1)
    m <- matrix(NA,nrow=nObs,ncol=2)
    colnames(m) <- c("cluster","child")
    
    m[,1] <- rep(1:nClust,each=nChildPerClust)
    m[,2] <- 1:nrow(m)
    m<-data.frame(cbind(m,x.cor))
    names(m) <- c("cluster","child",
                  'tx','x1','x2','x3','x4')
    
    y_t <- b0t + b1*m$x1[m$tx==1] + b2*m$x2[m$tx==1] +b3*m$x3[m$tx==1] + b4*m$x4[m$tx==1] + rnorm(nT,0,sigma.e)
    y_c <- b0c + b1*m$x1[m$tx==0] + b2*m$x2[m$tx==0] +b3*m$x3[m$tx==0] + b4*m$x4[m$tx==0] + rnorm(nC,0,sigma.e)
    
    y <- (c(y_t,y_c))
    
    m$y <- y + rep(rnorm(nClust,0,tau),each=nChildPerClust)
    
    ###############
    reg <- summary(lm(m$y~m$tx+m$x1+m$x2+m$x3+m$x4))
    Mtch[s,3] <- reg$coefficients[[2]]
    Mtch[s,4] <- reg$coefficients[[8]]
    
    I <- getICC(m$y,m$cluster)
    ICC[s,] <- c(I[[1]],I[[2]],I[[3]],I[[4]],I[[5]])
    
    ate.r <- summary(lm(m$y~m$tx))
    Mtch[s,1] <- ate.r$coefficients[[2]]
    Mtch[s,2] <- ate.r$coefficients[[4]]
    ###############PSM model
    Y <- m$y
    Tr <- m$tx
    glm1 <- glm(tx ~ x1 + x2 + x3 + x4, family = binomial, data = m)
    
    rr1 <- Match(Y = Y, Tr = Tr, X = glm1$fitted, estimand="ATE")
    rr2 <- Match(Y = Y, Tr = Tr, X = glm1$fitted, estimand="ATT")
    
    Mtch[s,5] <- rr1$est
    Mtch[s,6] <- rr1$se
    Mtch[s,7] <- rr2$est
    Mtch[s,8] <- rr2$se
    
    #IPW model
    m$pihat.log <- glm1$fitted
    #Calculate Weights
    
    m$ipw.ate <- ifelse(m$tx==1, 1/m$pihat.log,1/(1-m$pihat.log))
    #ATE Outcome Analysis
    
    design.ate <- svydesign(ids= ~1, weights= ~ipw.ate, dat=m)
    mod.ipw.ate <- summary(svyglm(y ~ tx, design=design.ate))
    
    Mtch[s,9] <- mod.ipw.ate$coefficients[[2]]#slope
    Mtch[s,10] <- mod.ipw.ate$coefficients[[4]]#sd
    #############################
    subcl <- matchit(tx ~ x1 + x2 + x3 +x4, data=m,
                     method='subclass',   subclass=5)
    #ATE Outcome Analysis
    sub.mtch<- summary(lm(y~ tx,weights=weights,data=match.data(subcl)))
    Mtch[s,11] <- sub.mtch$coefficients[[2]]#slope
    Mtch[s,12] <- sub.mtch$coefficients[[4]]#sd
    
    setTxtProgressBar(pb, s)
    
  }
  
  icc<- mean(ICC[,3])
  #################################naive
  naive.ate <- sum(Mtch[,1]/S)
  raw.bias.ate <- sum(Mtch[,1]/S)-g
  relative.bias.ate <- (((sum(Mtch[,1]/S)-g)/g) * 100)
  sd.ate <- sqrt(((S-1)^-1)*sum((Mtch[,1]-(sum(Mtch[,1]/S)))^2)) # more precise sd
  MSE.ate <- (sum(Mtch[,1]/S)-g)^2 + (sqrt(((S-1)^-1)*sum((Mtch[,1]-(sum(Mtch[,1]/S)))^2)))^2
  ###########ols
  ate.reg <- sum(Mtch[,3]/S)
  raw.bias <- sum(Mtch[,3]/S)-g
  relative.bias <- (((sum(Mtch[,3]/S)-g)/g) * 100)
  sd.reg <- sqrt(((S-1)^-1)*sum((Mtch[,3]-(sum(Mtch[,3]/S)))^2)) # more precise sd
  MSE.reg <- (sum(Mtch[,3]/S)-g)^2 + (sqrt(((S-1)^-1)*sum((Mtch[,3]-(sum(Mtch[,3]/S)))^2)))^2
  ####################psm
  psm <- sum(Mtch[,5]/S)
  raw.bias.psm <- sum(Mtch[,5]/S)-g
  relative.bias.psm <- (((sum(Mtch[,5]/S)-g)/g) * 100)
  sd.psm <- sqrt(((S-1)^-1)*sum((Mtch[,5]-(sum(Mtch[,5]/S)))^2))
  MSE.psm <- (sum(Mtch[,5]/S)-g)^2 + (sqrt(((S-1)^-1)*sum((Mtch[,5]-(sum(Mtch[,5]/S)))^2)))^2
  ########ipw
  ipw.ate <- sum(Mtch[,9]/S)
  raw.bias.ipw <- sum(Mtch[,9]/S)-g
  relative.bias.ipw <- (((sum(Mtch[,9]/S)-g)/g) * 100)
  sd.ipw <- sqrt(((S-1)^-1)*sum((Mtch[,9]-(sum(Mtch[,9]/S)))^2))
  MSE.ipw <- (sum(Mtch[,9]/S)-g)^2 + (sqrt(((S-1)^-1)*sum((Mtch[,9]-(sum(Mtch[,9]/S)))^2)))^2
  #####subclass
  sub <- sum(Mtch[,11]/S)
  raw.bias.sub <- sum(Mtch[,11]/S)-g
  relative.bias.sub <- (((sum(Mtch[,11]/S)-g)/g) * 100)
  sd.sub <- sqrt(((S-1)^-1)*sum((Mtch[,11]-(sum(Mtch[,11]/S)))^2))
  MSE.sub <- (sum(Mtch[,11]/S)-g)^2 + (sqrt(((S-1)^-1)*sum((Mtch[,11]-(sum(Mtch[,11]/S)))^2)))^2
  ###########psm.att
  psm.tt <- sum(Mtch[,7]/S)
  raw.bias.att <- sum(Mtch[,7]/S)-g
  relative.bias.att <- (((sum(Mtch[,7]/S)-g)/g) * 100)
  sd.att <- sqrt(((S-1)^-1)*sum((Mtch[,7]-(sum(Mtch[,7]/S)))^2))
  MSE.att <- (sum(Mtch[,7]/S)-g)^2 + (sqrt(((S-1)^-1)*sum((Mtch[,7]-(sum(Mtch[,7]/S)))^2)))^2
  
  result<-data.frame(cbind(naive.ate,raw.bias.ate,relative.bias.ate,sd.ate,MSE.ate,ate.reg,raw.bias,relative.bias,sd.reg,MSE.reg,
                           psm,raw.bias.psm,relative.bias.psm,sd.psm,MSE.psm,ipw.ate,
                           raw.bias.ipw,relative.bias.ipw,sd.ipw,MSE.ipw,sub,raw.bias.sub,relative.bias.sub,
                           sd.sub,MSE.sub,psm.tt,raw.bias.att,relative.bias.att,sd.att,MSE.att,icc))
  
  assign(paste0("result", tau), result)
  filename1 <- paste("result",tau,".csv", sep ="")
  write.csv( get(paste0("result",tau)) , file = filename1)
  
  assign(paste0("ICC", tau), ICC)
  filename2 <- paste("ICC",tau,".csv", sep ="")
  write.csv(get(paste0("ICC",tau)), file=filename2)
  
  assign(paste0("Mtch", tau), Mtch)
  filename3 <- paste("Mtch",tau,".csv", sep ="")
  write.csv(get(paste0("Mtch",tau)), file=filename3)
  rm(result)
  rm(ICC)
  rm(Mtch)
}
