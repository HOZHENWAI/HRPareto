source('rhrpar.r')
#Exact simulations
# Nrep=1000
# d=c(3)
# n=c(10,50,100,1000)
# alpha=c(0.5,1,1.2)
# RES=data.frame(n=rep(n,each=length(d)*length(alpha)),d=rep(d,each=length(alpha),times=length(n)),alpha=rep(alpha,length(d)*length(n)))
# for (i in 1:dim(RES)[1]){
#   print(i)
#   Q=diag(rep(1,RES$d[i]))-1/RES$d[i]
#   l=rep(-RES$alpha[i]/RES$d[i],RES$d[i])
#   y=replicate(Nrep,{
#     S=rhrpar(RES$n[i],Q,l)
#     hrparfit(S,init=list(Q=Q,l=l),form='concat')$Qlconcat
#   })
#   y=data.frame(t(y))
#   bias=apply(y,2,mean)-QltocQl(Q,l)
#   var=apply(y,2,var)
#   RES$bias_alpha[i]=-sum(cQltoQl(bias,RES$d[i])$l)
#   RES$var_alpha[i]=sum(cov(y)[(RES$d[i]*(RES$d[i]-1)/2+1):length(bias),(RES$d[i]*(RES$d[i]-1)/2+1):length(bias)])
#   RES$bias_Q11[i]=cQltoQl(bias,RES$d[i])$Q[1,1]
#   RES$var_Q11[i]=sum(cov(y)[1:(RES$d[i]-1),1:(RES$d[i]-1)])
# }

Nrep=1000
d=3
n=c(10,50,100,1000)
alpha=c(0.5,1,1.2)
RES=data.frame(n=rep(n,each=length(d)*length(alpha)),d=rep(d,each=length(alpha),times=length(n)),alpha=rep(alpha,length(d)*length(n)))
for (i in 1:dim(RES)[1]){
  print(i)
  Q=diag(rep(1,RES$d[i]))-1/RES$d[i]
  l=rep(-RES$alpha[i]/RES$d[i],RES$d[i])
  y=replicate(Nrep,{
    S=rhrpar(RES$n[i],Q,l)
    hrparfit_sym(S,init=c(Q[1,1],RES$alpha[i]))$estimate
  })
  y=data.frame(t(y))
  bias=apply(y,2,mean)-c(Q[1,1],RES$alpha[i])
  var=apply(y,2,var)
  RES$bias_alpha[i]=bias[2]
  RES$var_alpha[i]=var[2]
  RES$bias_Q11[i]=bias[1]
  RES$var_Q11[i]=var[1]
}
save.image("dim3sym_b.RData")