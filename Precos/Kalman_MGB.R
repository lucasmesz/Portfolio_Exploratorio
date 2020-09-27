library(readxl)
library(astsa)
library(zoo)
library(xts)
wti_df<- read_excel("IMPA_Mestrado_Profissional/Commodities e Futuros-2020/Exercicio_3/wti_impa.xlsx", 
                       col_types = c("date", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric"))
#View(wti)

#wti=xts(wti[,-1], order.by=as.Date(wti$Date, "%m/%d/%Y"))
#z= read.zoo(wti, header = TRUE, format = "%m/%d/%Y")

wti=xts(wti_df[,-1], order.by=as.Date(wti_df$Date))
###############
plot(wti$F2,type="l",lwd=2);
lines(wti$F9, type="l",col=2,lwd=2)
lines(wti$F14, type="l",col=3,lwd=2)
#############################

y = cbind(log(wti$F2),log(wti$F9),log(wti$F14));
num = nrow(y); 
input = rep(1,num)
A = array(rep(1,3), dim=c(3,1,num))
mu0 = -.35; Sigma0 = 1; Phi = 1
Temp = array(c(2/12,5/12,9/12), dim=c(3,1,1))

# Function to Calculate Likelihood
Linn = function(para){
  cQ = para[1]*sqrt(1/52) # sigma_w
  cR1 = para[2] # 11 element of chol(R)
  cR2 = para[3] # 22 element of chol(R)
  cR3 = para[4] # 33 element of chol(R)
  cR = matrix(c(cR1,0,0,0,cR2,0,0,0,cR3),3) # put the matrix together
  mi=para[5]
  #print(mi)
  c=para[6]
  Ups=(mi-c-0.5*cQ^2)*(1/52)
  r=para[7]
  Gam=(r-c)*Temp
  #print(Gam)
  #print(dim(y[1,]))
  kf = Kfilter1(num,y,A,mu0,Sigma0,Phi,Ups,Gam,cQ,cR,input)
  return(kf$like) }
# Estimation
init.par = c(.1,.1,.1,.1,.1,.1,.1) # initial values of parameters

(est = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE,
             control=list(trace=1,REPORT=1))) # output not shown




SE = sqrt(diag(solve(est$hessian)))
tstat=est$par/SE
# Display estimates
u = cbind(estimate=est$par, SE,tstat)
rownames(u)=c('sigw','cR1', 'cR2', 'cR3', 'mi','c','r')
colnames(u)=c('valor estimado','desvio-padrao','estatistica-t')
u

# Smooth (first set parameters to their final estimates)
cQ = est$par[1]
cR1 = est$par[2]
cR2 = est$par[3]
cR3 = est$par[4]
cR = matrix(c(cR1,0,0,0,cR2,0,0,0,cR3),3) # put the matrix together
(R = t(cR)%*%cR) # to view the estimated R matrix
mi=est$par[5]
c=est$par[6]
r=est$par[7]
Ups=(mi-c-0.5*cQ^2)*(1/52)
Gam=(r-c)*Temp


ks=Ksmooth1(num,y,A,mu0,Sigma0,Phi,Ups,Gam,cQ,cR,input)
#ks = Ksmooth1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
# Plot
xsm = ts(as.vector(exp(ks$xs)))
rmse = ts(sqrt(as.vector(ks$Ps)))
plot(xsm[0:200], col='black',ylab='Temperature Deviations',type='l')

lines(wti_df$F14[0:200], type="l",col='red',lwd=2)
lines(ts(as.vector(exp(ks$xp[0:200]))), type="l",col='blue',lwd=2)

xx = c(time(xsm), rev(time(xsm)))

yy = c(xsm-2*rmse, rev(xsm+2*rmse))
polygon(xx, yy, border=NA, col=gray(.6, alpha=.25))
lines(globtemp, type='o', pch=2, col=4, lty=6)
lines(globtempl, type='o', pch=3, col=3, lty=6)
