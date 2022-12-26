##########################
##########################

# Diagnóstico e envelope para o modelo Bernoulli
source("https://www.ime.unicamp.br/~cnaber/diag_Bern.r")
source("https://www.ime.unicamp.br/~cnaber/envel_Binom.r")


# Dados sobre germinaçãoo

library(xtable)
library(plotrix)
library(plyr)
library(gamlss)
library(qualityTools)

op <- options()

# sucesso
y <- c(10,23,23,26,17,5,53,55,32,46,10,8,10,8,23,0,3,22,15,32,3)

# total de replicações
m <- c(39,62,81,51,39,6,74,72,51,79,13,16,30,28,45,4,12,41,30,51,7)

# Dados
tabela <- cbind(c(y[1:5],0),c(m[1:5],0),y[6:11],m[6:11],c(y[12:16],0),c(m[12:16],0),c(y[17:21],0),c(m[17:21],0))
xtable(tabela)


# proporçãode sucessos
vp <- y/m
semente <- c(rep("O.aegyptiaca75",11),rep("O.aegyptiaca73",10))
semente <- as.factor(semente)
extrato <- c(rep("feijão",5),rep("pepino",6),rep("feijão",5),rep("pepino",5))
extrato <- as.factor(extrato)


# Análise descritiva (das proporções)

# por extrato
dpextra <- data.frame(vp,extrato)
rvp1 <- ddply(dpextra,.(extrato),summarise,media=mean(vp),dp=sqrt(var(vp)),vari=var(vp),cv=100*((sqrt(var(vp))/mean(vp))),n=length(vp))
xtable(rvp1)

# por semente
dpsem <- data.frame(vp,semente)
rvp2 <- ddply(dpsem,.(semente),summarise,media=mean(vp),dp=sqrt(var(vp)),vari=var(vp),cv=100*((sqrt(var(vp))/mean(vp))),n=length(vp))
xtable(rvp2)

# por semente x extração
dpextrasem <- data.frame(vp,extrato,semente)
rvp3 <- ddply(dpextrasem,.(extrato,semente),summarise,media=mean(vp),dp=sqrt(var(vp)),vari=var(vp),cv=100*((sqrt(var(vp))/mean(vp))),n=length(vp))
xtable(rvp3)

# por semente x extração (quantidade de replicas de Bernoulli)
dpextrasem <- data.frame(m,extrato,semente)
rvp4 <- ddply(dpextrasem,.(extrato,semente),summarise,soma=sum(m))

#options(digits=2)
#options(op)     # reset (all) initial options
#options("digits") 


# Tabelas com as proporções

tabelap <- rbind(c(rvp1$media[1],1-rvp1$media[1]),
                 c(rvp1$media[2],1-rvp1$media[2]),
                 c(rvp2$media[1],1-rvp2$media[1]),
                 c(rvp2$media[2],1-rvp2$media[2]))
xtable(tabelap)

rownames(tabelap) <- c("feijão","pepino","O.aegyptiaca73","O.aegyptiaca75")      


# Gráficos das proporções
plot(c(rvp1$media[1],1-rvp1$media[1],rvp1$media[2],1-rvp1$media[2]),axes=FALSE,ylab="proporções",xlab="taxa de germinação por tipo de extrato de raiz",cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19)
axis(1,c(1,2,3,4),labels=c("feijão-germinou","feijão-não germinou","pepino-germinou","pepino-não germinou"))
#axis(2,at=c(round(min(ta1),2),round(max(ta1),2)))
#axis(2,at=seq(round(min(ta1),2),round(max(ta1),2),2))
axis(2)


plot(c(rvp2$media[1],1-rvp2$media[1],rvp2$media[2],1-rvp2$media[2]),axes=FALSE,ylab="proporções",xlab="taxa de germinação por esp???cie",cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19)
axis(1,c(1,2,3,4),labels=c("O.aegyptiaca73-germinou","O.aegyptiaca73-não germinou","O.aegyptiaca75-germinou","O.aegyptiaca75-não germinou"))
#axis(2,at=c(round(min(ta1),2),round(max(ta1),2)))
#axis(2,at=seq(round(min(ta1),2),round(max(ta1),2),2))
axis(2)


# Gráfico de perfis
medias <- rvp3$media
mediasa <- medias
dp <- rvp3$dp
n <- rvp3$n
vn <- rvp4$soma
ez <- qnorm(0.975)

#plot(medias[1:2],axes=FALSE,ylim=c(0.15,0.80),cex.lab=1.5,xlab="tipo de semente",ylab="proporçãode sementes germinadas")
#axis(2,cex.axis=1.2)
#axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
#plotCI(medias[1:2],liw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[1:2]),uiw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[1:2]),pch=19,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=4,cex=1.2)
#lines(medias[1:2],lwd=2,col=4)
#plotCI(medias[3:4],liw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[3:4]),uiw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[3:4]),pch=23,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=2,pt.bg=2,cex=1.2)
#lines(medias[3:4],col=2,lwd=2)
#legend(1,0.35,col=c(4,2),lwd=c(2,2),pch=c(19,23),pt.bg=c(2,2),legend=c("feijão","pepino"),bty="n",cex=1.5)

plot(medias[1:2],axes=FALSE,ylim=c(0.00,0.80),cex.lab=1.5,xlab="tipo de semente",ylab="proporçãode sementes germinadas")
axis(2,cex.axis=1.2)
axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
plotCI(medias[1:2],liw=ez*dp[1:2]/sqrt(n[1:2]),uiw=ez*dp[1:2]/sqrt(n[1:2]),pch=19,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=4,cex=1.2)
lines(medias[1:2],lwd=2,col=4)
plotCI(medias[3:4],liw=ez*dp[3:4]/sqrt(n[3:4]),uiw=ez*dp[3:4]/sqrt(n[3:4]),pch=23,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=2,pt.bg=2,cex=1.2)
lines(medias[3:4],col=2,lwd=2)
legend(1.2,0.15,col=c(4,2),lwd=c(2,2),pch=c(19,23),pt.bg=c(2,2),legend=c("feijão","pepino"),bty="n",cex=1.5)



# modelo logito
model1 <- fit.model <- result <- glm(cbind(y,m-y)~semente+extrato+semente*extrato,family=binomial(link="logit"))
summary(result)
xtable(summary(result))
###  OBS > a iteração entre semente e extrato deu significativa, dessa forma, temos que considerar os dois fatores, mesmo que a varivel semente, individualmente, nap tenha dado significativa
###      > sementeO.aegyptiaca75:extratopepino = 0.7781: quando temos a especie 75 com o pepino, temos um aumento de 0.7781 na probabilidade de germinação
###      > sementeO.aegyptiaca75 = -0.1459: quando muda da especie 73 para a 75, temos uma queda de -0,15 na probabilidade de germinação
###      > extratopepino = 0.5401: quando muda do feijão para o pepino, temos um aumento de 0.5401 na probabilidade de germinação

desvio <- deviance(result)
pvdesv <- 1-pchisq(desvio,df=21-4)

ez <- qnorm(0.975)
pred <- predict(fit.model,type=c("response"),se.fit = TRUE) # valores preditos ; probabilidades estimadas e seu respectivo erro padrão
mupred <- unique(pred$fit) # valores preditos para cada classificação
mupredALL <- pred$fit # valores preditos para cada observação
sepredALL <- pred$se.fit # erro padrão para cada observação

mupred <- c(mupred[3],mupred[1],mupred[4],mupred[2])
semupred <- unique(pred$se.fit)
semupred <- c(semupred[3],semupred[1],semupred[4],semupred[2])

# IC para a media dos valores preditos de cada uma das classificações
liIC = unique(apply(cbind(0,mupred-ez*semupred),1,max))
lsIC = unique(apply(cbind(1,mupred+ez*semupred),1,min))

# IC para a media dos valores preditos de cada uma das observações
liICALL = apply(cbind(0,mupredALL-ez*sepredALL),1,max)
lsICALL = apply(cbind(1,mupredALL+ez*sepredALL),1,min)

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


# modelo probito
model2 <- fit.model <- result <- glm(cbind(y,m-y)~semente+extrato+semente*extrato,family=binomial(link="probit"))
summary(result)
xtable(summary(result))

diagBern(fit.model)
envelBinom(fit.model,"probit")

desvioLP <- deviance(result)
pvdesvLP <- 1-pchisq(desvioLP,df=21-4)

predLP <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLP <- unique(predLP$fit)
mupredALP <- pred$fit
sepredALP <- pred$se.fit

mupredLP <- c(mupredLP[3],mupredLP[1],mupredLP[4],mupredLP[2])
semupredLP <- unique(predLP$se.fit)
semupredLP <- c(semupredLP[3],semupredLP[1],semupredLP[4],semupredLP[2])

# IC para a media dos valores preditos de cada uma das classificações
liICLP = apply(cbind(0,mupredLP-ez*semupredLP),1,max)
lsICLP = apply(cbind(1,mupredLP+ez*semupredLP),1,min)

# IC para a media dos valores preditos de cada uma das observações
liICALP = (apply(cbind(0,mupredALP-ez*sepredALP),1,max))
lsICALP = (apply(cbind(1,mupredALP+ez*sepredALP),1,min))

AICLP <- AIC(fit.model)
BICLP <- BIC(fit.model)

rebeta <- (summary(fit.model))$coef
rebetaLP <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLP)
mrebeta <- rbind(mrebeta,rebetaLP)
covbetaLP <- vcov(fit.model)


# modelo cauchy
model3 <- fit.model <- result <- glm(cbind(y,m-y)~semente+extrato+semente*extrato,family=binomial(link="cauchit"))
summary(result)
xtable(summary(result))

diagBern(fit.model)
envelBinom(fit.model,"cauchit")

desvioLC<-deviance(result)
pvdesvLC <- 1-pchisq(desvioLC,df=21-4)


predLC <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLC <- unique(predLC$fit)
mupredALC <- pred$fit
sepredALC <- pred$se.fit

mupredLC <- c(mupredLC[3],mupredLC[1],mupredLC[4],mupredLC[2])
semupredLC <- unique(predLC$se.fit)
semupredLC <- c(semupredLC[3],semupredLC[1],semupredLC[4],semupredLC[2])

# IC para a media dos valores preditos de cada uma das classificações
liICLC = apply(cbind(0,mupredLC-ez*semupredLC),1,max)
lsICLC = apply(cbind(1,mupredLC+ez*semupredLC),1,min)

# IC para a media dos valores preditos de cada uma das observações
liICALC = (apply(cbind(0,mupredALC-ez*sepredALC),1,max))
lsICALC = (apply(cbind(1,mupredALC+ez*sepredALC),1,min))

AICLC <- AIC(fit.model)
BICLC <- BIC(fit.model)

rebeta <- (summary(fit.model))$coef
rebetaLC <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLC)
mrebeta <- rbind(mrebeta,rebetaLC)
covbetaLC <- vcov(fit.model)


# modelo cloglog
model4 <- fit.model <- result <- glm(cbind(y,m-y)~semente+extrato+semente*extrato,family=binomial(link="cloglog"))
summary(result)
xtable(summary(result))

desvioLCLL<- deviance(result)
pvdesvLCLL <- 1-pchisq(desvioLCLL,df=21-4)

diagBern(fit.model)
envelBinom(fit.model,"cloglog")

predLCLL <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLCLL <- unique(predLCLL$fit)
mupredALCLL <- pred$fit
sepredALCLL <- pred$se.fit

mupredLCLL <- c(mupredLCLL[3],mupredLCLL[1],mupredLCLL[4],mupredLCLL[2])
semupredLCLL <- unique(predLCLL$se.fit)
semupredLCLL <- c(semupredLCLL[3],semupredLCLL[1],semupredLCLL[4],semupredLCLL[2])

# IC para a media dos valores preditos de cada uma das classificações
liICLCLL = apply(cbind(0,mupredLCLL-ez*semupredLCLL),1,max)
lsICLCLL = apply(cbind(1,mupredLCLL+ez*semupredLCLL),1,min)

# IC para a media dos valores preditos de cada uma das observações
liICALCLL = (apply(cbind(0,mupredALCLL-ez*sepredALCLL),1,max))
lsICALCLL = (apply(cbind(1,mupredALCLL+ez*sepredALCLL),1,min))

AICLCLL <- AIC(fit.model)
BICLCLL <- BIC(fit.model)

rebeta <- (summary(fit.model))$coef
rebetaLCLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLCLL)
mrebeta <- rbind(mrebeta,rebetaLCLL)
xtable(mrebeta)
covbetaLCLL <- vcov(fit.model)


# AIC e BIC
AICBIC <- rbind(cbind(AICLL,AICLP,AICLC,AICLCLL),cbind(BICLL,BICLP,BICLC,BICLCLL))
colnames(AICBIC) <- c("logito","probito","cauchito","cloglog")
rownames(AICBIC) <- c("AIC","BIC")
xtable(t(AICBIC))


# Desvio absolutos médios
damLL <- mean(abs(y/m-mupredALL))
damLP <- mean(abs(y/m-mupredALP))
damLC <- mean(abs(y/m-mupredALC))
damLCLL <- mean(abs(y/m-mupredLCLL))

AICBICDAM <- rbind(cbind(AICLL,AICLP,AICLC,AICLCLL),cbind(BICLL,BICLP,BICLC,BICLCLL),cbind(damLL,damLP,damLC,damLCLL),cbind(desvio,desvioLP,desvioLC,desvioLCLL),cbind(pvdesv,pvdesvLP,pvdesvLC,pvdesvLCLL))
colnames(AICBICDAM) <- c("logito","probito","cauchito","cloglog")
rownames(AICBICDAM) <- c("AIC","BIC","DAM","Desvio","p-valor desvio")
xtable(t(AICBICDAM))


# Envelopes dos quatro modelos
par(mfrow=c(2,2))

envelBinom(model1,"logit")
title("ligação logito")

envelBinom(model2,"probit")
title("ligação probito")

envelBinom(model3,"cauchit")
title("ligação cauchito")

envelBinom(model4,"cloglog")
title("ligação cloglog")


# Valores preditos (proporções médias) para os quatro modelos
par(mfrow=c(2,2))

plot(medias[1:2],axes=FALSE,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
lines(medias[3:4],col=1,type="p",cex=1.2,pch=17)
plotCI(mupred[1:2],li=liIC[1:2],ui=lsIC[1:2],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=19,col=2,add=TRUE)
lines(mupred[1:2],col=2,type="l",cex=1.2,pch=17,lwd=2)
plotCI(mupred[3:4],li=liIC[3:4],ui=lsIC[3:4],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=17,col=2,add=TRUE)
lines(mupred[3:4],col=2,type="l",cex=1.2,pch=17,lwd=2)
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.9,1.2,col=c(1,1),lty=c(0,0),lwd=c(0,0),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (observada)","pepino (observada)"),bty="n",cex=1.2)
legend(1.5,1.2,col=c(2,2),lty=c(1,1),lwd=c(2,2),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
title("ligação logito")


plot(medias[1:2],axes=FALSE,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
lines(medias[3:4],col=1,type="p",cex=1.2,pch=17)
plotCI(mupredLP[1:2],li=liICLP[1:2],ui=lsICLP[1:2],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=19,col=2,add=TRUE)
lines(mupredLP[1:2],col=2,type="l",cex=1.2,pch=17,lwd=2)
plotCI(mupredLP[3:4],li=liICLP[3:4],ui=lsICLP[3:4],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=17,col=2,add=TRUE)
lines(mupredLP[3:4],col=2,type="l",cex=1.2,pch=17,lwd=2)
title("ligação probito")
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.9,1.2,col=c(1,1),lty=c(0,0),lwd=c(0,0),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (observada)","pepino (observada)"),bty="n",cex=1.2)
legend(1.5,1.2,col=c(2,2),lty=c(1,1),lwd=c(2,2),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)


plot(medias[1:2],axes=FALSE,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
lines(medias[3:4],col=1,type="p",cex=1.2,pch=17)
plotCI(mupredLC[1:2],li=liICLC[1:2],ui=lsICLC[1:2],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=19,col=2,add=TRUE)
lines(mupredLC[1:2],col=2,type="l",cex=1.2,pch=17,lwd=2)
plotCI(mupredLC[3:4],li=liICLC[3:4],ui=lsICLC[3:4],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=17,col=2,add=TRUE)
lines(mupredLC[3:4],col=2,type="l",cex=1.2,pch=17,lwd=2)
title("ligação cauchito")
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.9,1.2,col=c(1,1),lty=c(0,0),lwd=c(0,0),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (observada)","pepino (observada)"),bty="n",cex=1.2)
legend(1.5,1.2,col=c(2,2),lty=c(1,1),lwd=c(2,2),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)


plot(medias[1:2],axes=FALSE,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
axis(2,cex.axis=1.2)
axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
lines(medias[3:4],col=1,type="p",cex=1.2,pch=17)
plotCI(mupredLCLL[1:2],li=liICLCLL[1:2],ui=lsICLCLL[1:2],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=19,col=2,add=TRUE)
lines(mupredLCLL[1:2],col=2,type="l",cex=1.2,pch=17,lwd=2)
plotCI(mupredLCLL[3:4],li=liICLCLL[3:4],ui=lsICLCLL[3:4],cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="ln(dose)",ylab="proporçãode insetos mortos",pch=17,col=2,add=TRUE)
lines(mupredLCLL[3:4],col=2,type="l",cex=1.2,pch=17,lwd=2)
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.9,1.2,col=c(1,1),lty=c(0,0),lwd=c(0,0),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (observada)","pepino (observada)"),bty="n",cex=1.2)
legend(1.5,1.2,col=c(2,2),lty=c(1,1),lwd=c(2,2),pch=c(19,17),pt.bg=c(2,2),legend=c("feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
title("ligação cloglog")


# Valores preditos (proporções individuais) para os quatro modelos
par(mfrow=c(2,2))
plot(vp,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
#axis(2,cex.axis=1.2)
#axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
plotCI(mupredALL,li=liICALL,ui=lsICALL,pch=19,col=2,add=TRUE)
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.05,1.3,col=c(1,2),lty=c(0,0),lwd=c(0,0),pch=c(19,19),pt.bg=c(2,2),legend=c("observada","predita"),bty="n",cex=1.2)
title("ligação logit")


plot(vp,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
#axis(2,cex.axis=1.2)
#axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
plotCI(mupredALP,li=liICALP,ui=lsICALP,pch=17,col=2,add=TRUE)
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.05,1.3,col=c(1,2),lty=c(0,0),lwd=c(0,0),pch=c(19,19),pt.bg=c(2,2),legend=c("observada","predita"),bty="n",cex=1.2)
title("ligação probit")


plot(vp,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
#axis(2,cex.axis=1.2)
#axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
plotCI(mupredALC,li=liICALC,ui=lsICALC,pch=17,col=2,add=TRUE)
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.05,1.3,col=c(1,2),lty=c(0,0),lwd=c(0,0),pch=c(19,19),pt.bg=c(2,2),legend=c("observada","predita"),bty="n",cex=1.2)
title("ligação cauchito")


plot(vp,ylim=c(0.2,1.2),cex.lab=1.2,xlab="tipo de semente",ylab="proporçãode sementes germinadas",pch=19,cex=1.2,col=1)
#axis(2,cex.axis=1.2)
#axis(1,1:2,c("O.aegyptiaca73","O.aegyptiaca75"),cex.axis=1.2)
plotCI(mupredALCLL,li=liICALCLL,ui=lsICALCLL,pch=17,col=2,add=TRUE)
#legend(1,1.2,col=c(1,1,2,2),lty=c(0,0,1,1),lwd=c(0,0,2,2),pch=c(19,17,19,17),pt.bg=c(2,2,2,2),legend=c("feijão (observada)","pepino (observada)","feijão (preditiva)","pepino (preditiva)"),bty="n",cex=1.2)
legend(0.05,1.3,col=c(1,2),lty=c(0,0),lwd=c(0,0),pch=c(19,19),pt.bg=c(2,2),legend=c("observada","predita"),bty="n",cex=1.2)
title("ligação cloglog")