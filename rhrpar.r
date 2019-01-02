if(!exists("pmvnorm",mode="function")) library('mvtnorm')
library('stats')
#library('numDeriv')
rhrpar=function(n,Q,l,a=rep(1,length(l))){
  #compute alpha the exponent of the Pareto distribution
  alpha=-sum(l)
  #extract dimension
  d=length(l);
  #compute the corresponding parameter for a HRPareto with margin equal to 1_d
  l=l-as.numeric(Q%*%log(a))
  #Generate n observations of a real Pareto distributed random variable with parameter alpha
  R=runif(n)^(-1/alpha)
  #Compute p the vector of probabilities for the componentwise maxima,i.e. p[i] is the probability that the i-th component is maximal
  p=NULL
  for (i in 1:d){
    p[i]=(det(solve(as.matrix(Q[-i,-i]))))^(1/2)*exp((1/2)*as.numeric(t(l[-i])%*%solve(Q[-i,-i])%*%l[-i]))*pmvnorm(upper=rep(0,d-1),algorithm=Miwa(),mean=as.numeric(solve(Q[-i,-i])%*%l[-i]),sigma=solve(Q[-i,-i]))
  }
  p=p/sum(p)
  #Recursive sampling procedure
  y=NULL;
  #Generate the indices of maximum for each observation
  jm=sample(1:d,n,replace=TRUE,prob=p);
  for (i in 1:n){
    #Initialise G as a d dimentional vector of zeros
    g=rep(1,d-1)
    #Recursively generate random gaussian variables
    J=(1:d)[-jm[i]]
    Sig=solve(Q[-jm[i],-jm[i]])
    mu=Sig%*%l[-jm[i]]
    while(max(g)>0){
      g=rmvnorm(1,mu,Sig)
    }
    G=rep(0,d)
    G[(1:d)[-jm[i]]]=g
    #Compute Theta
    theta=exp(G)
    #Generate the sample following a H?ssler-Reiss with margin 1_d
    Z=R[i]*theta  
    #Use transformation properties to generate the final sample with margin a
    Z=a*Z
    #Add sequentially each samples
    y=rbind(y,as.numeric(Z))
  }
  #Return
  return(y)
}



Ca=function(Q,l,a=rep(1,length(l))){
  y=(2*pi)^((length(l)-1)/2)/sum(-l)
  s=0;
  for (i in 1:length(l)){
    s=s+(a[i])^(sum(l))*det(solve(Q[-i,-i]))^(1/2)*exp(0.5*as.numeric(l[-i]%*%solve(Q[-i,-i])%*%l[-i]))*pmvnorm(upper=log(a[-i]/a[i]),mean=as.numeric(solve(Q[-i,-i])%*%l[-i]),sigma=solve(Q[-i,-i]),algorithm=Miwa())
  }
  return(as.numeric(y*s))
}


statshrpar=function(S){
  d=dim(S)[2]
  n=dim(S)[1]
  St=list(V=matrix(rep(0,d^2),ncol=d),M=rep(0,d))
  m=log(S)%*%rep(1,d)/d
  m=matrix(rep(m,times=d,each=1),ncol=d)
  St$V=-1/(2*n)*t(log(S)-m)%*%(log(S)-m)
  St$M=(1/n)*apply(log(S),2,sum)
  return(St)
}

cQltoQl=function(cQl,d){
  Q=matrix(0,d,d)
  l=cQl[(d*(d-1)/2+1):length(cQl)]
  Q[upper.tri(Q)]=cQl[1:(d*(d-1)/2)]
  Q=Q+t(Q)
  diag(Q)=apply(Q,1,sum)
  diag(Q)=-diag(Q)
  return(list(Q=Q,l=l))
}

QltocQl=function(Q,l){
  return(c(Q[upper.tri(Q)],l))
}


nhrparlogl=function(cQl,S,d,a=rep(1,d),concatQl=TRUE,concatS=TRUE){
  #If cQl is concatened so that the cQl is a vector where the d(d-1)/2 first component are the upper diag of Q and the rest are the components of l
  if (concatQl){
    parameters=cQltoQl(cQl,d)
  }else{ #other cQl must be a list with cQl$Q and cQl$l
    parameters=cQl
  }
  #Compute the statistics for easier computing if needed
  if(!concatS){#The default format for S the sufficient statistics
    S=statshrpar(S)
  }
  if ((sum(eigen(parameters$Q)$values>(10^(-3)))!=(d-1)) | sum(parameters$l)>0) {return(10^3)}
  else {
    K=-(sum(S$V*parameters$Q))-sum(S$M*parameters$l)+log(Ca(parameters$Q,parameters$l,a))
    if (is.na(K) | is.infinite(K)){
      print(Ca(parameters$Q,parameters$l,a))
    }
    return(K)
  }
}

wrapper_sym=function(par,d,a,S){
  #par[1]=Q[1,1]
  #par[2]=alpha
  #S is the statistic
  c=par[1]
  alpha=par[2]
  Q=diag(rep(c+c/(d-1),d))-c/(d-1)
  l=rep(-alpha/d,d)
  return(nhrparlogl(list(Q=Q,l=l),S,d,a=a,concatQl=FALSE))
}

hrparfit=function(S,...,a=rep(1,dim(S)[2]),form="Ql"){
  # retrieve optional argument
  op=list(...)
  #n the size the sample and d the dimension
  n=dim(S)[1]
  d=dim(S)[2]
  #if some initialisation is available, either init has to be a (d(d-1)/2+d) vector or a list that has to have two component, one named Q and the other named l 
  if (is.null(op$init)){
    tol=0
    init=-runif(min=10^-2,d*(d-1)/2+d)
    I=cQltoQl(init,d)
    while (nhrparlogl(init,S,d,a,concatS=FALSE)==10^3){
      init=-runif(d*(d-1)/2+d)
      I=cQltoQl(init,d)
    }
  }else{
    init=op$init
    if(exists("init",mode="list")){
      init=QltocQl(init$Q,init$l)
    }
  }
  #Compute the exhaustive statistics
  st=statshrpar(S)
  #Since the problem is strictly convex, the nlm routine is enough to compute the minimum for the negative loglikelihood
  y=nlm(nhrparlogl,init,S=st,d=d,a=a,hessian=FALSE,ndigit=6,steptol=1e-10)
  #y=optim(init,nhrparlogl,method="Nelder-Mead",S=st,d=d,a=a)
  #convert the result into desired format
  x=NULL
  if (form=="concat"){x$Qlconcat=y$estimate}
  if (form=="Ql"){x=cQltoQl(y$estimate,d)}
  x$minimum=y$minimum
  #x$gradient=y$gradient
  #x$hessian=y$hessian
  #computing information criterion
  x$AIC=2*(d-1+d*(d-1)/2)+2*y$minimum
  x$BIC=log(n)*(d-1+d*(d-1)/2)+2*y$minimum
  return(x)
  #x return the hessian for the upper diagonal component of Q and l
  #more work has to be done to get confidence interval for the diagonal of Q
}

hrparfit_sym=function(S,...,a=rep(1,dim(S)[2]),form="Ql"){
  op=list(...)
  #n the size the sample and d the dimension
  n=dim(S)[1]
  d=dim(S)[2]
  st=statshrpar(S)
  #init must be specified
  init=op$init
  y=nlm(wrapper_sym,init,S=st,d=d,a=a,hessian=FALSE,ndigit=6,steptol=1e-10)
  return(y)
}
##############################################Generalized HR

rghrpar=function(n,alpha,Q,l,a=rep(1,length(l))){
  #general a H?ssler-Reiss pareto sample of length l
  Y=rhrpar(n,Q,l-Q%*%diag(alpha)%*%log(a))
  #The transformation properties are used to convert the HR sample into a generalised HRPareto sample with exponent alpha in R^d
  Y=a*(Y^t(matrix(rep((-sum(l-Q%*%diag(alpha)%*%log(a))/alpha),n),nrow=length(l))))
  #return the sample, each rows is a realisation of a random variable with distribution HRPar(alpha,Q,l) and margin a
  return(Y)
}

caQltoaQl=function(caQl,d){
  res=NULL
  Q=matrix(0,d,d)
  Q[upper.tri(Q)]=caQl[(d+1):(d+d*(d-1)/2)]
  Q=Q+t(Q)
  diag(Q)=apply(Q,1,sum)
  diag(Q)=-diag(Q) #kerQ=vect(1d)
  res$Q=Q
  res$l=c(-1-sum(caQl[(d+1+d*(d-1)/2):length(caQl)]),caQl[(d+1+d*(d-1)/2):length(caQl)])#t(l)1d=-1
  res$alpha=caQl[1:d] #alpha in (0,infty)^d
  return(res)
}

aQltocaQl=function(alpha,Q,l){
  caQl=NULL
  d=length(l)
  caQl=rep(0,d+d-1+d*(d-1)/2)
  caQl[1:d]=alpha
  caQl[(d+1):(d+d*(d-1)/2)]=Q[upper.tri(Q)]
  caQl[(d+1+d*(d-1)/2):length(caQl)]=l[2:length(l)]
  return(caQl)
}

cgQltoQl=function(cgQl,d){
  res=NULL
  Q=matrix(0,d,d)
  Q[upper.tri(Q)]=cgQl[1:(d*(d-1)/2)]
  Q=Q+t(Q)
  diag(Q)=apply(Q,1,sum)
  diag(Q)=-diag(Q) #kerQ=vect(1d)
  res$Q=Q
  res$l=c(-1-sum(cgQl[(d*(d-1)/2+1):length(cgQl)]),cgQl[(d*(d-1)/2+1):length(cgQl)])
  return(res)
}

QltocgQl=function(Q,l){
  d=length(l)
  cgQl=rep(0,d*(d-1)/2+d-1)
  cgQl[1:(d*(d-1)/2)]=Q[upper.tri(Q)]
  cgQl[(d*(d-1)/2+1):length(cgQl)]=l[2:length(l)]
  return(cgQl)
}

gstatshrpar=function(S){
  return(list(V=-t(log(S))%*%log(S)/(2*dim(S)[1]),M=apply(log(S),2,mean))) 
}

nghrparlogl=function(pcaQl,S,d,a=rep(1,d),concatQl=TRUE,concatS=TRUE){
  #is concatQl is TRUE, pcaQl must be converted into a list with 3 component: alpha, Q ,l
  if (concatQl){
    parameters=caQltoaQl(pcaQl,d)}
  else{
    parameters=pcaQl
  }
  #if concatS=FALSE, there is a need to compute the sufficient statistics
  if (!concatS){S=gstatshrpar(S)}
  #compute normalizing term
  if (sum(eigen(parameters$Q)$values>(10^(-3)))!=(d-1)| sum(parameters$alpha>0)!=d){return(10^3)}
  else {
    C=prod(parameters$alpha)^(-1)*Ca(parameters$Q,parameters$l,a^(parameters$alpha))
    lC=log(C)
    #Compute the scalar product <T(z),h(\theta)>
    D_a=diag(parameters$alpha)
    Q=parameters$Q
    l=parameters$l
    r=sum(S$V*(D_a%*%Q%*%D_a))+t(l)%*%D_a%*%S$M
    #return the negative log likelihood
    res=lC-r
    return(res)}
}

split_nghrparlogl1=function(alpha,cgQl,S,d,a=rep(1,d)){                   #this looks really slow though 
  par=cgQltoQl(cgQl,d)
  return(nghrparlogl(pcaQl = aQltocaQl(alpha,par$Q,par$l),S,d,a))
}

split_nghrparlogl2=function(cgQl,alpha,S,d,a=rep(1,d)){                   #too many conversions
  par=cgQltoQl(cgQl,d)
  return(nghrparlogl(pcaQl = aQltocaQl(alpha,par$Q,par$l),S,d,a))
}

############################alternative minisation
alternative=function(ialpha,iQ,il,St,tol=10^-2,stepmax=100,d=length(ialpha),a=rep(1,d)){ 
  res=NULL
  step=1
  stop=100
  #initialisation
  halpha=ialpha
  hgQl=QltocgQl(iQ,il)
  while (step<stepmax & stop>tol){
    T1=nlm(f=split_nghrparlogl1,halpha,cgQl=hgQl,S=St,d=d,a=a,ndigit=6)
    halpha=T1$estimate
    T2=nlm(f=split_nghrparlogl2,hgQl,alpha=halpha,S=St,d=d,a=a,ndigit=6)
    hgQl=T2$estimate
    step=step+1
    stop=abs(T1$minimum-T2$minimum)
  }
  hQtl=cgQltoQl(hgQl,d)
  res$estimate=aQltocaQl(halpha,hQtl$Q,hQtl$l)
  res$minimum=T2$minimum
  return(res)
}

ghrparfit=function(S,a=rep(1,dim(S)[2]),method='nlm',form='Ql'){
  #initialisation + retrieving dimension and sample length
  res=NULL
  d=dim(S)[2]
  n=dim(S)[1]
  ######################compute stats
  
  St=gstatshrpar(S)
  
  ####################initialising using moment estimators
  
  N=apply(S>t(matrix(rep(a,n),nrow=d)),2,mean)   #N is the proportion of excedance in each dimension
  O=apply((S>t(matrix(rep(a,n),nrow=d)))*log(S),2,mean) #O is the mean of the excedance in each dimension
  hatalpha=(O/N-log(a))^(-1)  #since the marginal follow pareto (a_i,alpha_i), taking the ratio of the probability of excedance over the expectation of the log excedance returns 1/alpha_i
  initQl=list(Q=diag(rep(1,d))-1/d,l=rep(-sum(hatalpha)/d,d))
  hatQl=cgQltoQl(nlm(split_nghrparlogl2,QltocgQl(initQl$Q,initQl$l),alpha=hatalpha,S=St,d=d,a=a,ndigit = 6)$estimate,d) #computing argmax_Ql
  
  #################### Once the parameter are initialised, to obtain maximum likelihood estimator, one can use the strict concavity of the loglikehood function around the true value by using optimisation routine for the general function
  if (method=='nlm'){
    pcaQl=aQltocaQl(hatalpha,hatQl$Q,hatQl$l)
    ans=nlm(nghrparlogl,pcaQl,S=St,d=d,a=a,hessian=TRUE,ndigit=6)
  }
  
  ####################alternatively, one can use an alternating minimisation scheme, the gradient and hessian has to be computed manually
  
  if(method=='alternative'){
    ans=alternative(hatalpha,hatQl$Q,hatQl$l,St,a=a)
    ans$gradient=grad(nghrparlogl,ans$estimate,S=St,d=d)
    ans$hessian=hessian(nghrparlogl,ans$estimate,S=St,d=d)
  }
  
  #########################formatting results depending on the arguments, by default return a list
  if (form=='Ql'){
    res=caQltoaQl(ans$estimate,d)
    res$minimum=ans$minimum
    res$gradient=ans$gradient
    res$hessian=ans$hessian
  }else{
    if (form=='concat'){
      res=ans
      names(res)[2]='aQlconcat'
    }
  }
  
  #########################computing information criteria
  res$AIC=2*(d-1+d*(d-1)/2)+2*res$minimum
  res$BIC=log(n)*(d-1+d*(d-1)/2)+2*res$minimum
  #########################return
  return(res)
}

hrpartest=function(S,a=rep(1,dim(S)[2]),alpha_level=0.05){
  d=dim(S)[2]
  test=NULL
  L_0=hrparfit(S=S,a=a,form="concat")$minimum
  L_1=ghrparfit(S=S,a=a,form="concat")$minimum
  Delta=L_0-L_1
  test$threshold_right=qchisq(1-alpha_level,df=d-1)
  test$statistics=2*Delta
  test$p_value=1-pchisq(test$statistics,df=d-1)
  test$df=d-1
  test$hyp='Null hypothesis H_0: alpha_1=alpha_2=...=alpha_d'
  return(test)
}

