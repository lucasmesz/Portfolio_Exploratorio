
library(astsa)
plot(globtemp,type="l",ylim=c(-.5,.8),lwd=2);
lines(globtempl, type="l",col=2,lwd=2)
y = cbind(globtemp,globtempl); num = nrow(y); input = rep(1,num)
A = array(rep(1,2), dim=c(2,1,num))
mu0 = -.35; Sigma0 = 1; Phi = 1
# Function to Calculate Likelihood
Linn = function(para){
  cQ = para[1] # sigma_w
  cR1 = para[2] # 11 element of chol(R)
  cR2 = para[3] # 22 element of chol(R)
  cR12 = para[4] # 12 element of chol(R)
  cR = matrix(c(cR1,0,cR12,cR2),2) # put the matrix together
  drift = para[5]
  kf = Kfilter1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
  return(kf$like) }
# Estimation
init.par = c(.1,.1,.1,0,.05) # initial values of parameters
(est = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE,
             control=list(trace=1,REPORT=1))) # output not shown




SE = sqrt(diag(solve(est$hessian)))
tstat=est$par/SE
# Display estimates
u = cbind(estimate=est$par, SE,tstat)
rownames(u)=c('sigw','cR11', 'cR22', 'cR12', 'drift')
colnames(u)=c('valor estimado','desvio-padrao','estatistica-t')
u

# Smooth (first set parameters to their final estimates)
cQ = est$par[1]
cR1 = est$par[2]
cR2 = est$par[3]
cR12 = est$par[4]
cR = matrix(c(cR1,0,cR12,cR2), 2)
(R = t(cR)%*%cR) # to view the estimated R matrix
drift = est$par[5]
ks = Ksmooth1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
# Plot
xsm = ts(as.vector(ks$xs), start=1880)
rmse = ts(sqrt(as.vector(ks$Ps)), start=1880)
plot(xsm, ylim=c(-.6, 1), ylab='Temperature Deviations')
xx = c(time(xsm), rev(time(xsm)))
yy = c(xsm-2*rmse, rev(xsm+2*rmse))
polygon(xx, yy, border=NA, col=gray(.6, alpha=.25))
lines(globtemp, type='o', pch=2, col=4, lty=6)
lines(globtempl, type='o', pch=3, col=3, lty=6)
