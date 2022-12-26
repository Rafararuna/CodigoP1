######################################
######################################
######################################
#          Regressão Binária
######################################
######################################
######################################

library(xtable)
library(plotrix)
library(plyr)
op <- options()

# Diagnóstico e envelope para o modelo Bernoulli
source("https://www.ime.unicamp.br/~cnaber/diag_Bern.r")
source("https://www.ime.unicamp.br/~cnaber/envel_Bern.r")
source("https://www.ime.unicamp.br/~cnaber/envel_Binom.r")



######################################
######################################
######################################
# Exemplo da mortalidade dos besouros
######################################
######################################


x = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839) # as doses aplicadas
m = c(59, 60, 62, 56, 63, 59, 62, 60) # numero de insetos expostos as doses aplicadas em x
y = c(6, 13, 18, 28, 52, 53, 61, 60) # numero de insetos expostos que morreram, dos m que foram expostos, as doses aplicadas em x
proba = y/m
plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19)


#result<- glm(y/m~x,family=binomial(link="logit"),weights=m)


xc <- x-mean(x) # x corrigido

fit.model <- result <- glm(cbind(y,m-y)~xc,family=binomial(link="logit"))
summary(result)
xtable(summary(result))

desvio <- deviance(result)
pvdesv <- 1-pchisq(desvio,df=6) # > 0,05, então não há evidências contra o modelo

diagBern(fit.model)
envelBinom(fit.model,"logit")

AICLL <- AIC(fit.model)
BICLL <- BIC(fit.model)

ez <- qnorm(0.975)
rebeta <- (summary(result))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)
mrebeta <- rbind(rebetaLL)
covbetaLL <- vcov(fit.model)


# Estimativas das proporções com IC's (contas à mão)

## ez <- qnorm(0.975)
## vbeta <- as.numeric(coef(result))
## m.X <- model.matrix(result)
## m.X <- matrix(as.numeric(m.X),8,2)
## p0<- fitted(result)
## aux1 <- p0*(1-p0)
## covbeta <- vcov(result)
## covbeta <-  matrix(as.numeric(covbeta),2,2)
## epbeta <- sqrt(diag(covbeta))
## mPsi <- aux1^2*(m.X%*%covbeta%*%t(m.X))
## epprob <- sqrt(diag(mPsi))
## mICprob <- cbind(p0-ez*epprob,p0+ez*epprob) 

## plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=1,ylim=c(0,1))
## plotCI(x,p0,li=mICprob[,1],ui=mICprob[,2],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,add=TRUE)
## legend(1.70,0.6,legend=c("observada","predita"),pch=c(1,19),bty="n",cex=1.2)


# Via função predict

ez <- qnorm(0.975)
pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit # probabilidades prediras ; valor de Mi
semupred <- pred$se.fit # erros padrão preditos de Mi
liIC = apply(cbind(0,mupred-ez*semupred),1,max)
lsIC = apply(cbind(1,mupred+ez*semupred),1,min)
plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=1,ylim=c(0,1))
plotCI(x,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,add=TRUE)
legend(1.70,0.9,legend=c("observada","predita"),pch=c(1,19),bty="n",cex=1.2)



# Outras funções de ligação

## PROBITO

fit.model <- result <- glm(cbind(y,m-y)~xc,family=binomial(link="probit"))
summary(result)
xtable(summary(result))

desvioLP <- deviance(result)
pvdesvLP <- 1-pchisq(desvioLP,df=6)

diagBern(fit.model)
envelBinom(fit.model,"probit")

predLP <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLP <- predLP$fit
semupredLP <- predLP$se.fit
liICLP = apply(cbind(0,mupredLP-ez*semupredLP),1,max)
lsICLP = apply(cbind(1,mupredLP+ez*semupredLP),1,min)

AICLP <- AIC(fit.model)
BICLP <- BIC(fit.model)

rebeta <- (summary(fit.model))$coef
rebetaLP <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLP)
mrebeta <- rbind(mrebeta,rebetaLP)
covbetaLP <- vcov(fit.model)

## CAUCHY

fit.model <- result <- glm(cbind(y,m-y)~xc,family=binomial(link="cauchit"))
summary(result)
xtable(summary(result))

desvioLC <- deviance(result)
pvdesvLC <- 1-pchisq(desvioLC,df=6)

diagBern(fit.model)
envelBinom(fit.model,"cauchit")

predLC <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLC <- predLC$fit
semupredLC <- predLC$se.fit
liICLC = apply(cbind(0,mupredLC-ez*semupredLC),1,max)
lsICLC = apply(cbind(1,mupredLC+ez*semupredLC),1,min)

AICLC <- AIC(fit.model)
BICLC <- BIC(fit.model)

rebeta <- (summary(fit.model))$coef
rebetaLC <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLC)
mrebeta <- rbind(mrebeta,rebetaLC)
covbetaLC <- vcov(fit.model)

## COMPLEMENTO LOG-LOG

fit.model<-result <- glm(cbind(y,m-y)~xc,family=binomial(link="cloglog"))
summary(result)
xtable(summary(result))

desvioLCLL <- deviance(result)
pvdesvLCLL <- 1-pchisq(desvioLCLL,df=6)

diagBern(fit.model)
envelBinom(fit.model,"cloglog")

predLCLL <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLCLL <- predLCLL$fit
semupredLCLL <- predLCLL$se.fit
liICLCLL = apply(cbind(0,mupredLCLL-ez*semupredLCLL),1,max)
lsICLCLL = apply(cbind(1,mupredLCLL+ez*semupredLCLL),1,min)

AICLCLL <- AIC(fit.model)
BICLCLL <- BIC(fit.model)

rebeta <- (summary(fit.model))$coef
rebetaLCLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLCLL)
mrebeta <- rbind(mrebeta,rebetaLCLL)

xtable(mrebeta)
covbetaLCLL <- vcov(fit.model)



#
AICBIC<- rbind(cbind(AICLL,AICLP,AICLC,AICLCLL),cbind(BICLL,BICLP,BICLC,BICLCLL))
colnames(AICBIC)<- c("logito","probito","cauchito","cloglog")
rownames(AICBIC)<- c("AIC","BIC")
xtable(t(AICBIC))



# Desvio absolutos médios

damLL <- mean(abs(y/m-mupred))
damLP <- mean(abs(y/m-mupredLP))
damLC <- mean(abs(y/m-mupredLC))
damLCLL <- mean(abs(y/m-mupredLCLL))

AICBICDAM<- rbind(cbind(AICLL,AICLP,AICLC,AICLCLL),cbind(BICLL,BICLP,BICLC,BICLCLL),cbind(damLL,damLP,damLC,damLCLL),cbind(desvio,desvioLP,desvioLC,desvioLCLL),cbind(pvdesv,pvdesvLP,pvdesvLC,pvdesvLCLL))
colnames(AICBICDAM)<- c("logito","probito","cauchito","cloglog")
rownames(AICBICDAM)<- c("AIC","BIC","DAM","Desvio","p-valor desvio")
xtable(t(AICBICDAM))



# Predição de todos os modelos

plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,ylim=c(0,1))
plotCI(x,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col=2)
plotCI(x,mupredLP,li=liICLP,ui=lsICLP,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col=3)
plotCI(x,mupredLC,li=liICLC,ui=lsICLC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col=4)
plotCI(x,mupredLCLL,li=liICLCLL,ui=lsICLCLL,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col="yellow3")
legend(1.70,0.8,legend=c("observado","predito - logito","predito - probito","predito - cauchito","predito - clolog"),pch=c(19,17,17,17,17),col=c(1,2,3,4,"yellow3"),bty="n",cex=1.2)

par(mfrow=c(2,2))
plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,ylim=c(0,1),main="modelo logito")
plotCI(x,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,col=2,add=TRUE)
legend(1.68,0.8,legend=c("observado","predito"),pch=c(19,17),col=c(1,2),bty="n",cex=1.2)
plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,ylim=c(0,1),main="modelo probito")
plotCI(x,mupredLP,li=liICLP,ui=lsICLP,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,col=3,add=TRUE)
legend(1.68,0.8,legend=c("observado","predito"),pch=c(19,17),col=c(1,3),bty="n",cex=1.2)
plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,ylim=c(0,1),main="modelo cauchito")
plotCI(x,mupredLC,li=liICLC,ui=lsICLC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,col=4,add=TRUE)
legend(1.68,0.8,legend=c("observado","predito"),pch=c(19,17),col=c(1,4),bty="n",cex=1.2)
plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,ylim=c(0,1),main="modelo cloglog")
plotCI(x,mupredLCLL,li=liICLCLL,ui=lsICLCLL,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,col="yellow3",add=TRUE)
legend(1.68,0.8,legend=c("observado","predito"),pch=c(19,17),col=c(1,"yellow3"),bty="n",cex=1.2)


#
#plot(x,proba,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,ylim=c(0,1),xlim=c(1.68,1.95))
#plotCI(x+0.006,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=19,add=TRUE,col=2)
#plotCI(x+0.012,mupredLP,li=liICLP,ui=lsICLP,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col=3)
#plotCI(x+0.018,mupredLC,li=liICLC,ui=lsICLC,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col=4)
#plotCI(x+0.024,mupredLCLL,li=liICLCLL,ui=lsICLCLL,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporção de insetos mortos",pch=17,add=TRUE,col=5)



# Estimativas de doses letais

pdose <- c(0.5,0.7,0.8,0.99)


# DS

DSLL <- c(log(pdose/(1-pdose))-rebetaLL[1,1])/rebetaLL[2,1] + mean(x)

DSLP <- c(qnorm(pdose)-rebetaLP[1,1])/rebetaLP[2,1] + mean(x)

#DSLC <- c(qt(pdose,df=1)-rebetaLC[1,1])/rebetaLC[2,1] + mean(x)
DSLC <- c(qcauchy(pdose)-rebetaLC[1,1])/rebetaLC[2,1] + mean(x)

DSLCLL <- c(log(-log(1-pdose))-rebetaLCLL[1,1])/rebetaLCLL[2,1] + mean(x)


# EP

EPDSLL <- sqrt(diag(cbind(-1/rebetaLL[2,1],cbind((1/rebetaLL[2,1]^2)*(rebetaLL[1,1]-log(pdose/(1-pdose)))))%*%covbetaLL%*%t(cbind(-1/rebetaLL[2,1],cbind((1/rebetaLL[2,1]^2)*(rebetaLL[1,1]-log(pdose/(1-pdose))))))))
EPDSLP <- sqrt(diag(cbind(-1/rebetaLP[2,1],cbind((1/rebetaLP[2,1]^2)*(rebetaLP[1,1]-qnorm(pdose))))%*%covbetaLP%*%t(cbind(-1/rebetaLP[2,1],cbind((1/rebetaLP[2,1]^2)*(rebetaLP[1,1]-qnorm(pdose)))))))
EPDSLC <- sqrt(diag(cbind(-1/rebetaLC[2,1],cbind((1/rebetaLC[2,1]^2)*(rebetaLC[1,1]-qcauchy(pdose))))%*%covbetaLC%*%t(cbind(-1/rebetaLC[2,1],cbind((1/rebetaLC[2,1]^2)*(rebetaLC[1,1]-qcauchy(pdose)))))))
EPDSLCLL <- sqrt(diag(cbind(-1/rebetaLCLL[2,1],cbind((1/rebetaLCLL[2,1]^2)*(rebetaLCLL[1,1]-log(-log(1-pdose)))))%*%covbetaLCLL%*%t(cbind(-1/rebetaLCLL[2,1],cbind((1/rebetaLCLL[2,1]^2)*(rebetaLCLL[1,1]-log(-log(1-pdose))))))))


# IC

ICDSLL <- c(DSLL-ez*EPDSLL,DSLL+ez*EPDSLL)
ICDSLP <- c(DSLP-ez*EPDSLP,DSLP+ez*EPDSLP)
ICDSLC <- c(DSLC-ez*EPDSLC,DSLC+ez*EPDSLC)
ICDSLCLL <- c(DSLCLL-ez*EPDSLCLL,DSLCLL+ez*EPDSLCLL)

#
aux <- length(pdose)
par(mfrow=c(1,1))
plotCI(pdose,DSLL,li=ICDSLL[1:aux],ui=ICDSLL[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=1,xlim=c(0.5,1.05),ylim=c(1.75,2))
plotCI(pdose+0.01,DSLP,li=ICDSLP[1:aux],ui=ICDSLP[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=2,add=TRUE)
plotCI(pdose+0.02,DSLC,li=ICDSLC[1:aux],ui=ICDSLC[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=3,add=TRUE)
plotCI(pdose+0.03,DSLCLL,li=ICDSLCLL[1:aux],ui=ICDSLCLL[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=4,add=TRUE)
legend(0.55,2,legend=c("logito","probito","cauchito","clolog"),pch=c(17,17,17,17),col=c(1,2,3,4),bty="n",cex=1.2)

#
plotCI(pdose,DSLL,li=ICDSLL[1:aux],ui=ICDSLL[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=1,xlim=c(0.5,1.05),ylim=c(1.75,2.7))
plotCI(pdose+0.01,DSLP,li=ICDSLP[1:aux],ui=ICDSLP[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=2,add=TRUE)
plotCI(pdose+0.02,DSLC,li=ICDSLC[1:aux],ui=ICDSLC[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=3,add=TRUE)
plotCI(pdose+0.03,DSLCLL,li=ICDSLCLL[1:aux],ui=ICDSLCLL[(aux+1):(2*aux)],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="proporção de interesse",ylab="dose necessária",pch=17,col=4,add=TRUE)
legend(0.55,2.5,legend=c("logito","probito","cauchito","clolog"),pch=c(17,17,17,17),col=c(1,2,3,4),bty="n",cex=1.2)