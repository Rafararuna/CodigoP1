######################################
######################################
######################################
#          Regress?o Bin?ria
######################################
######################################
######################################

library(xtable)
library(plotrix)
library(plyr)
op <- options()

# Diagn?stico e envelope para o modelo Bernoulli:
source("https://www.ime.unicamp.br/~cnaber/diag_Bern.r")
source("https://www.ime.unicamp.br/~cnaber/envel_Bern.r")



###########################################
###########################################
#         Dados sobre vasoconstri??o
###########################################
###########################################

m.dados <- scan("https://www.ime.unicamp.br/~cnaber/pregibon.dat",what=list(0,0,0))
v.y <- cbind(m.dados[[1]]) # vasocontri??o
v.vol <- cbind(m.dados[[2]]) # volume
v.raz <- cbind(m.dados[[3]])  # raz?o
v.lvol <- log(v.vol) # log do volume
v.lraz <- log(v.raz) # log da raz?o
n <- length(v.y)


# Dados
mdados < -cbind(v.y,v.vol,v.raz)
xtable(mdados)


# Gr?ficos de dispers?o

par(mfrow=c(2,2))
plot(v.vol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="volume",ylab="vasocontri??????o",pch=19)
plot(v.raz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="raz???o",ylab="vasocontri??????o",pch=19)
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??????o",pch=19)
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz???o)",ylab="vasocontri??????o",pch=19)

par(mfrow=c(1,1))
plot(v.vol,v.lraz,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="ln(raz???o)",pch=19)
text(v.vol,v.lraz,labels=v.y,pos=4,cex=1.2)


# An?lise descritiva do log das covari?veis

## Transformando os fatores segundo parametriza??o casela de refer?ncia
yfac<- factor(v.y,levels=c("0","1"))#

## data.frame para usar a fun??o ddply

dvc<-data.frame(v.lvol,v.lraz,yfac)
rvol<-ddply(dvc,.(yfac),summarise,media=mean(v.lvol),median=quantile(v.lvol,0.5),dp=sqrt(var(v.lvol)),vari=var(v.lvol),cv=100*((sqrt(var(v.lvol))/mean(v.lvol))),minimo=min(v.lvol),maximo=max(v.lvol))
rraz<-ddply(dvc,.(yfac),summarise,media=mean(v.lraz),median=quantile(v.lraz,0.5),dp=sqrt(var(v.lraz)),vari=var(v.lraz),cv=100*((sqrt(var(v.lraz))/mean(v.lraz))),minimo=min(v.lraz),maximo=max(v.lraz))
mres <- data.frame(t(rvol),t(rraz))
options(digits=2)
amres <-data.frame(t(rvol),t(rraz))
options(op)     # reset (all) initial options
options("digits")
xtable(amres)


# Box plots
par(mfrow=c(1,2))
boxplot(v.vol~v.y,cex=1.2,cex.axis=1.2,cex.label=1.2,xlab="ocorr?ncia de vasoconstri?ao",ylab="ln(volume)",names=c("n?o","sim"))
boxplot(v.lraz~v.y,cex=1.2,cex.axis=1.2,cex.label=1.2,xlab="ocorr?ncia de vasoconstri??o",ylab="ln(raz?o)",names=c("n?o","sim"))


# An?lise inferencial

result <- fit.model<-glm(v.y~v.lvol+v.lraz,family=binomial("logit"))
summary(fit.model)

rebeta <- (summary(result))$coef

diagBern(fit.model) # diagnostico do modelo
envelBern(fit.model,"logit")

ez <- qnorm(0.975)

rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)
AICLL <- AIC(fit.model)
BICLL <- BIC(fit.model)
mrebeta <- rbind(rebetaLL)


# M?dias preditas

pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit # probabilidades/m?dias preditas
semupred <- pred$se.fit # erros padr?o preditos das probabilidades preditas
liIC = apply(cbind(0,mupred-ez*semupred),1,max) # limite inferior
lsIC = apply(cbind(1,mupred+ez*semupred),1,min) # limite superios
#lines(v.lvol,mupred,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??????o",pch=19,type="p",col=2)

par(mfrow=c(1,2))
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??o",pch=19,xlim=c(-1.3,1.31))
plotCI(v.lvol,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-1.49,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz?o)",ylab="vasocontri??o",pch=19)
plotCI(v.lraz,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-3,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)


# Valor predito para a vari?vel resposta

ypred <- rbinom(n,1,mupred)
par(mfrow=c(1,2))
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??o",pch=19,main="ocorr?ncias de vasocontri??o observadas e preditas pelo modelo")
lines(v.lvol,ypred,type="p",cex=1.2,pch=17,col=2)
legend(0,0.6,legend=c("observado","predito"),pch=c(19,17),col=c(1,2),bty="n",cex=1.2)
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz?o)",ylab="vasocontri??o",pch=19,main="ocorr?ncias de vasocontri??o observadas e preditas pelo modelo")
lines(v.lraz,ypred,type="p",cex=1.2,pch=17,col=2)
legend(-1,0.6,legend=c("observado","predito"),pch=c(19,17),col=c(1,2),bty="n",cex=1.2)



# Como obter intervalos de confian?a para as raz?es de chances? 
# 1) Fazer um IC para o par?metro original e depois calcular o IC para a tranforma??o
# 2) M?todo Delta
# 3) Simula??o/Reamostragem

## 1) 

### Obten??o dos IC's assint?ticas para as fun??es de interesse

vbeta <- as.numeric(coef(result))
covbeta <- vcov(result)
covbeta <-  matrix(as.numeric(covbeta),3,3)
epbeta <- sqrt(diag(covbeta))

### Transforma??o

mICbeta <- cbind(vbeta-1.96*epbeta,vbeta+1.96*epbeta) # limites inferiores e superiores dos betas
mIC1 <- rbind(cbind(exp(mICbeta[1,1])/(1+exp(mICbeta[1,1])),exp(mICbeta[1,2])/(1+exp(mICbeta[1,2]))),
              cbind(exp(mICbeta[2,1]),exp(mICbeta[2,2])),cbind(exp(mICbeta[3,1]),exp(mICbeta[3,2]))) # intervalos da raz?o de chances para cada beta
#OBS: o problema da transforma??o ? que, depois que a gente aplica ela, n?o tem como mais garantir um n?vel de signific?ncia de 95% ainda


## 2) 

### M?todo delta

#### fun??es de interesse
g1 <- exp(vbeta[1])/(1+exp(vbeta[1])) 
g2 <- exp(vbeta[2]) 
g3 <- exp(vbeta[3]) 

edev <- deviance(result)

#### gradientes
psi1 <- cbind(exp(vbeta[1])/((1+exp(vbeta[1]))^2),0,0) # primeiro gradiente
psi2 <- cbind(0,exp(vbeta[2]),0) # segundo gradiente
psi3 <- cbind(0,0,exp(vbeta[3])) # terceiro gradiente

mPsi <- rbind(psi1,psi2,psi3) # media do psi
epfbeta <- sqrt(diag(mPsi%*%covbeta%*%t(mPsi))) # erro padr?o

mIC2 <- rbind(cbind(g1-1.96*epfbeta[1],g1+1.96*epfbeta[1]),
              cbind(g2-1.96*epfbeta[2],g2+1.96*epfbeta[2]),cbind(g3-1.96*epfbeta[3],g3+1.96*epfbeta[3]))
# OBS: os betas est?o nos Reais, ent?o a exponencial dos betas deviam estar entre 0 e 1 ;
#      por?m, percebe-se a presen?a de valores negativos, isso acontece porque assumimos normalidade para os betas, entao eles estao de - inf at? + inf


## 3)

### Simula??o

#### Obten??o num?rica da distribui??o da raz?o de chances

veta <- vbeta[1]+v.lvol*vbeta[2]+v.lraz*vbeta[3] # preditor linear dos parametros
vmu <- exp(veta)/(1+exp(veta)) # probabilidade estimada
ys <- rbinom(n,1,vmu) # valor predito
results <- glm(ys~v.lvol+v.lraz,family=binomial(link="logit")) # constroi um novo modelo
vbetas <- coef(results) # acha os novos valores para betas
mestrc <- cbind(exp(vbetas[1])/(1+exp(vbetas[1])),exp(vbetas[2]),exp(vbetas[3]))
for(j in 2:500)
{ # calculo das fun??es de interesse
  ys <- rbinom(n,1,vmu)
  results<- glm(ys~v.lvol+v.lraz,family=binomial(link="logit"))
  vbetas <- coef(results)
  mestrc <- rbind(mestrc,cbind(exp(vbetas[1])/(1+exp(vbetas[1])),exp(vbetas[2]),exp(vbetas[3])))
}

mIC3 <- cbind(apply(mestrc,2,quantile,0.025),apply(mestrc,2,quantile,0.975))


# xtable(cbind(mIC1,mIC2,mIC3))
est.tau <- rbind(exp(vbeta[1])/(1+exp(vbeta[1])),exp(vbeta[2]),exp(vbeta[3]))
xtable(cbind(est.tau ,mIC1,mIC2,mIC3))



######################################################################################



# Outras fun??es de liga??o

## Probito

fit.model <- glm(v.y~v.lvol+v.lraz,family=binomial("probit"))
diagBern(fit.model) # semelhante ao logito
envelBern(fit.model,"probit")
AICLP <- AIC(fit.model)
BICLP <- BIC(fit.model)
rebeta <- (summary(fit.model))$coef
rebetaLP <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLP)
mrebeta <- rbind(mrebeta,rebetaLP)

# M?dias preditas

predLP <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLP <- predLP$fit
semupredLP <- predLP$se.fit
liICLP = apply(cbind(0,mupredLP-ez*semupredLP),1,max)
lsICLP = apply(cbind(1,mupredLP+ez*semupredLP),1,min)

## Cauchito

fit.model <- glm(v.y~v.lvol+v.lraz,family=binomial("cauchit"))
diagBern(fit.model)
envelBern(fit.model,"cauchit")
AICLC <- AIC(fit.model)
BICLC <- BIC(fit.model)
rebeta <- (summary(fit.model))$coef
rebetaLC <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLC)
mrebeta <- rbind(mrebeta,rebetaLC)

# M?dias preditas

predLC <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLC <- predLC$fit
semupredLC <- predLC$se.fit
liICLC = apply(cbind(0,mupredLC-ez*semupredLC),1,max)
lsICLC = apply(cbind(1,mupredLC+ez*semupredLC),1,min)

## Cloglog

fit.model <- glm(v.y~v.lvol+v.lraz,family=binomial("cloglog"))
diagBern(fit.model)
envelBern(fit.model,"cloglog")
AICLCLL <- AIC(fit.model)
BICLCLL <- BIC(fit.model)
rebeta <- (summary(fit.model))$coef
rebetaLCLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLCLL)
mrebeta <- rbind(mrebeta,rebetaLCLL)

xtable(mrebeta)

# M?dias preditas

predLCLL <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLCLL <- predLCLL$fit
semupredLCLL <- predLCLL$se.fit
liICLCLL = apply(cbind(0,mupredLCLL-ez*semupredLCLL),1,max)
lsICLCLL = apply(cbind(1,mupredLCLL+ez*semupredLCLL),1,min)


# AIC e BIC

AICBIC <- rbind(cbind(AICLL,AICLP,AICLC,AICLL),cbind(BICLL,BICLP,BICLC,BICLL))
colnames(AICBIC) <- c("logito","probito","cauchito","cloglog")
rownames(AICBIC) <- c("AIC","BIC")
xtable(t(AICBIC))

# OBS: nao podemos comparar apenas olhando esses resultados, temos que fazer a compara??o grafica tbm


#######################################



# Poder preditivo de todos os modelos

par(mfrow=c(2,2))
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??o",pch=19,xlim=c(-1.3,1.31),main="liga??o logito")
plotCI(v.lvol,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-1.49,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??o",pch=19,xlim=c(-1.3,1.31),main="liga??o probito")
plotCI(v.lvol,mupredLP,li=liICLP,ui=lsICLP,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-1.49,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??o",pch=19,xlim=c(-1.3,1.31),main="liga??o cauchito")
plotCI(v.lvol,mupredLC,li=liICLC,ui=lsICLC,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-1.49,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lvol,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(volume)",ylab="vasocontri??o",pch=19,xlim=c(-1.3,1.31),main="liga??o complemento log-log")
plotCI(v.lvol,mupredLCLL,li=liICLCLL,ui=lsICLCLL,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-1.49,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)

par(mfrow=c(2,2))
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz?o)",ylab="vasocontri??o",pch=19,main="liga??o logito")
plotCI(v.lraz,mupred,li=liIC,ui=lsIC,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-3,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz?o)",ylab="vasocontri??o",pch=19,main="liga??o probito")
plotCI(v.lraz,mupredLP,li=liICLP,ui=lsICLP,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-3,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz?o)",ylab="vasocontri??o",pch=19,main="liga??o cauchito")
plotCI(v.lraz,mupredLC,li=liICLC,ui=lsICLC,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-3,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)
plot(v.lraz,v.y,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(raz?o)",ylab="vasocontri??o",pch=19,main="liga??o complemento log-log")
plotCI(v.lraz,mupredLCLL,li=liICLCLL,ui=lsICLCLL,cex=1.2,cex.axis=1.2,cex.lab=1.2,col=2,add=TRUE,pch=19)
legend(-3,0.9,legend=c("observado","predita"),pch=c(19,19),col=c(1,2),bty="n",cex=1.2)

# Desvio absolutos m?dios

damLL <- mean(abs(v.y-mupred))
damLP <- mean(abs(v.y-mupredLP))
damLC <- mean(abs(v.y-mupredLC))
damLCLL <- mean(abs(v.y-mupredLCLL))

AICBICDAM <- rbind(cbind(AICLL,AICLP,AICLC,AICLCLL),cbind(BICLL,BICLP,BICLC,BICLCLL),cbind(damLL,damLP,damLC,damLCLL))
                                                                                          colnames(AICBICDAM)<- c("logito","probito","cauchito","cloglog")
                                                                                          rownames(AICBICDAM)<- c("AIC","BIC","DAM")
                                                                                          xtable(t(AICBICDAM))


