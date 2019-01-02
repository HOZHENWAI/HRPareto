#Scatter plot

source('rhrpar.r')
#options globales
densite=function(x1,y1,c,l1,l2){
   z=c(x1,y1)
   if(x1<1 & y1<1)
   {
     return(0)
   }
   else
   {
   Q=matrix(c(c,-c,-c,c),nrow=2)
   l=c(l1,l2)
   val=exp(-0.5 *log(z)%*%Q%*%log(z)+l%*%log(z))/(x1*y1*Ca(Q,l))
   return(val)
   }
}
#png('rplot.png') # sauvegarde automatique sur rplot.ps - on peut remplacer postscript/ps par png pour avoir un ficher png
n=200  #taille echantillon
par(mfrow=c(2,3)) #nombre de graphes
d=2 #dimension

#premier plot  : 
c=1
alpha=1
delta=0

Q=matrix(c(c,-c,-c,c),nrow=2)
l=c(-alpha/2 +delta, -alpha/2 -delta)
S1=rhrpar(n,Q,l)
#deuxième plot :
c=1/2
alpha=1
delta=0

Q=matrix(c(c,-c,-c,c),nrow=2)
l=c(-alpha/2 +delta, -alpha/2 -delta)
S2=rhrpar(n,Q,l)
#troisième plot :
c=1
alpha=1
delta=alpha/2

Q=matrix(c(c,-c,-c,c),nrow=2)
l=c(-alpha/2 +delta, -alpha/2 -delta)
S3=rhrpar(n,Q,l)
#quatrième plot
c=1
alpha=2
delta=0

Q=matrix(c(c,-c,-c,c),nrow=2)
l=c(-alpha/2 +delta, -alpha/2 -delta)
S4=rhrpar(100,Q,l)
#cinquième :
c=1/2
alpha=2
delta=0

Q=matrix(c(c,-c,-c,c),nrow=2)
l=c(-alpha/2 +delta, -alpha/2 -delta)
S5=rhrpar(n,Q,l)
#sixième plot :
c=1
alpha=2
delta=0.5

Q=matrix(c(c,-c,-c,c),nrow=2)
l=c(-alpha/2 +delta, -alpha/2 -delta)
S6=rhrpar(n,Q,l)

#affichage
xmax=max(S1[,1],S2[,1],S3[,1],S4[,1],S5[,1],S6[,1])
xmin=min(S1[,1],S2[,1],S3[,1],S4[,1],S5[,1],S6[,1])
ymax=max(S1[,2],S2[,2],S3[,2],S4[,2],S5[,2],S6[,2])
ymin=min(S1[,2],S2[,2],S3[,2],S4[,2],S5[,2],S6[,2])
xlim=c(round(xmin,digits=3),xmax)
ylim=c(round(ymin,digits=3),ymax)


x=c(seq(0.05,1,by=0.05),seq(2,10,by=1),seq(20,100,by=10),seq(200,min(1000,xmax),by=100),seq(min(xmax,2000),xmax,by=1000))
y=c(seq(0.05,1,by=0.05),seq(2,10,by=1),seq(20,100,by=10),seq(200,min(1000,ymax),by=100),seq(min(ymax,2000),ymax,by=1000))

library(latex2exp)
z1=outer(x,y,Vectorize(densite),c=1,l1=-0.5,l2=-0.5)#a changer manuellement
plot(S1,xlab="", ylab="",xlim=xlim,ylim=ylim,log='xy', main=TeX('$\\alpha = 1,\ \\beta = 1,\ \\gamma = 0$'))
#image(x,y,z1,add=TRUE)

z2=outer(x,y,Vectorize(densite),c=0.5,l1=-0.5,l2=-0.5)
plot(S2,xlab="", ylab="",xlim=xlim,ylim=ylim,log='xy', main=TeX('$\\alpha = 1,\ \\beta = 0.5,\ \\gamma = 0$'))
#image(x,y,z2,add=TRUE)

z3=outer(x,y,Vectorize(densite),c=1,l1=0,l2=-1)
plot(S3,xlab="", ylab="",xlim=xlim,ylim=ylim,log='xy', main=TeX('$\\alpha = 1,\ \\beta = 1,\ \\gamma = 0.5$'))
#image(x,y,z3,add=TRUE)

z4=outer(x,y,Vectorize(densite),c=1,l1=-1,l2=-1)
plot(S4,xlab="", ylab="",xlim=xlim,ylim=ylim,log='xy',main=TeX('$\\alpha = 2,\ \\beta = 1,\ \\gamma = 0$'))
#image(x,y,z4,add=TRUE)

z5=outer(x,y,Vectorize(densite),c=0.5,l1=-1,l2=-1)
plot(S5,xlab="", ylab="",xlim=xlim,ylim=ylim,log='xy',main=TeX('$\\alpha = 2,\ \\beta = 0.5,\ \\gamma = 0$'))
#image(x,y,z5,add=TRUE)

z6=outer(x,y,Vectorize(densite),c=1,l1=0,l2=-2)
plot(S6,xlab="", ylab="",xlim=xlim,ylim=ylim,log='xy', main=TeX('$\\alpha = 2,\ \\beta = 1,\ \\gamma = 0.5$'))
#image(x,y,z6,add=TRUE)



#dev.off()
#second affichage

pdf('rplot_densite.pdf') # sauvegarde automatique sur rplot.ps - on peut remplacer postscript/ps par png pour avoir un ficher png
#taille echantillon
par(mfrow=c(2,3))

image(x,y,z1,xlim=c(0.05,xmax),ylim=c(0.05,ymax),log='xy')
mtext("c=1, delta= 0 , alpha= 1")

image(x,y,z2,xlim=c(0.05,xmax),ylim=c(0.05,ymax),log='xy')
mtext("c=1/2, delta= 0 , alpha= 1")

image(x,y,z3,xlim=c(0.05,xmax),ylim=c(0.05,ymax),log='xy')
mtext("c=1, delta= 1/2 , alpha= 1")

image(x,y,z4,xlim=c(0.05,xmax),ylim=c(0.05,ymax),log='xy')
mtext("c=1, delta= 0 , alpha= 2")

image(x,y,z5,xlim=c(0.05,xmax),ylim=c(0.05,ymax),log='xy')
mtext("c=1/2, delta= 0 , alpha= 2")

image(x,y,z6,xlim=c(0.05,xmax),ylim=c(0.05,ymax),log='xy')
mtext("c=1, delta= 1 , alpha= 2")
dev.off()

#extrapolation densit� par noyau
#premier plot  : 
# library(MASS)
# c=1
# alpha=1
# delta=0
# 
# Q=matrix(c(c,-c,-c,c),nrow=2)
# l=c(-alpha/2 +delta, -alpha/2 -delta)
# S1=rhrpar(10000,Q,l)
# d1=kde2d(S1[,1],S1[,2])
# image(d1,log='xy')
# #persp(d1)
# #deuxième plot :
# c=1/2
# alpha=1
# delta=0
# 
# Q=matrix(c(c,-c,-c,c),nrow=2)
# l=c(-alpha/2 +delta, -alpha/2 -delta)
# S2=rhrpar(10000,Q,l)
# #troisième plot :
# c=1
# alpha=1
# delta=alpha/2
# 
# Q=matrix(c(c,-c,-c,c),nrow=2)
# l=c(-alpha/2 +delta, -alpha/2 -delta)
# S3=rhrpar(10000,Q,l)
# #quatrième plot
# c=1
# alpha=2
# delta=0
# 
# Q=matrix(c(c,-c,-c,c),nrow=2)
# l=c(-alpha/2 +delta, -alpha/2 -delta)
# S4=rhrpar(10000,Q,l)
# #cinquième :
# c=1/2
# alpha=2
# delta=0
# 
# Q=matrix(c(c,-c,-c,c),nrow=2)
# l=c(-alpha/2 +delta, -alpha/2 -delta)
# S5=rhrpar(10000,Q,l)
# #sixième plot :
# c=1
# alpha=2
# delta=alpha/2
# 
# Q=matrix(c(c,-c,-c,c),nrow=2)
# l=c(-alpha/2 +delta, -alpha/2 -delta)
# S6=rhrpar(10000,Q,l)
