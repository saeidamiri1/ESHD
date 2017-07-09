sortx<-function(x){
  q0<-quantile(x,probs=c(.1,0.9))
  return(x[x>q0[1]&x<q0[2]])
}


mean1<-function(x) mean(sortx(x))
quantile1<-function(x,p) quantile(sortx(x),probs=p)

GeneSub<-function(Xtrai,Ytrai,B){
  resi<-resi2<- resi0<-NULL
  sran0<-m.lasso0 <- m.lasso<-list()
  i1<-1
  i2<-0
  end0<-FALSE
  while(i1<=B){
    i2<-i2+1
    sran0[[i1]]<-sample(dim(Xtrai)[2],length(Ytrai)-2,replace=FALSE)
    XX<-data.frame(Ytrai,Xtrai[,sran0[[i1]]])
    mol0<-lm(Ytrai~., XX)
    aa<-step(mol0,trace = FALSE,direction="backward")
    fval<-summary(aa)$fstatistic
    if(i1/i2<.20) {print("It seems the data is very sparse") ; end0<-TRUE; break}
    if(1-pf(fval[1],fval[2],fval[3])>0.2) next
    #    m.lasso[[i]] <- glmnet(as.matrix(Xtrai[,sran0[[i]]]),Ytrai,alpha=alpha0)
    i1<-i1+1
  }
  if(end0==TRUE) res<- NA else res<-sran0
  return(res) 
}

GeneSub2<-function(Xtrai,Ytrai,B){
  resi<-resi2<- resi0<-NULL
  sran0<-m.lasso0 <- m.lasso<-list()
  i1<-1
  i2<-0
  end0<-FALSE
  while(i1<=B){
    i2<-i2+1
    sran0[[i1]]<-sample(dim(Xtrai)[2],length(Ytrai),replace=FALSE)
    XX<-data.frame(Ytrai,Xtrai[,sran0[[i1]]])
    # mol0<-lm(Ytrai~., XX)
    #aa<-step(mol0,trace = FALSE,direction="backward")
    #fval<-summary(aa)$fstatistic
    # if(i1/i2<.20) {print("It seems the data is very sparse") ; end0<-TRUE; break}
    #  if(1-pf(fval[1],fval[2],fval[3])>0.2) next
    #    m.lasso[[i]] <- glmnet(as.matrix(Xtrai[,sran0[[i]]]),Ytrai,alpha=alpha0)
    i1<-i1+1
  }
  if(end0==TRUE) res<- NA else res<-sran0
  return(res) 
}




MODELS<-function(Xtrai,Ytrai,B,alpha0,sran00,Ncore){
  cl<-makeCluster(Ncore)
  registerDoParallel(cl) 
  clusterExport(cl, varlist = c
                ("Xtrai","Ytrai","glmnet","B","alpha0","sran00"), envir=environment())
  m.lasso<-parLapply(cl, 1:B, function(exponent)   glmnet(as.matrix(Xtrai[,sran00[[exponent]]]),Ytrai,alpha=alpha0))
  stopCluster(cl)
  return(list(sran00,m.lasso))
}



rslasso<-function(Xtrai,Ytrai,B=20, alpha0=1,Ncore=2,Fasy=TRUE ){
  
  if(Fasy=TRUE)   sran00<-GeneSub(Xtrai,Ytrai,B)
  else sran00<-GeneSub2(Xtrai,Ytrai,B)
  MO0<-MODELS(Xtrai,Ytrai,B,alpha0,sran00,Ncore)
  return(MO0)
}


prslasso<-function(MO0,Xtest){
  sran0<-MO0[[1]]
  m.lasso<-MO0[[2]]
  
  resii0<-NULL
  for(i in 1:(length(MO0[[1]]))){
    y.pred.lasso <- predict(m.lasso[[i]],as.matrix(Xtest[,sran0[[i]]]))
    resii0<-cbind(apply(y.pred.lasso,1,mean1)  ,resii0)
  }
  out0<-cbind(apply(resii0,1,mean1),apply(resii0,1,quantile1,.025),apply(resii0,1,quantile1,.975))
  colnames(out0)<-c('Estimate',"lower CI","upper CI")
  return(out0)
}
