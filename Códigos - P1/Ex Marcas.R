###########################################
###########################################
#         Dados sobre preferência de a
#         automóveis
###########################################
###########################################


library(xtable)
library(plotrix)
library(plyr)
op <- options()


# Diagnósstico e envelope para o modelo Bernoulli
source("https://www.ime.unicamp.br/~cnaber/diag_Bern.r")
source("https://www.ime.unicamp.br/~cnaber/envel_Bern.r")

# Funções de Análise do Desvio e do teste CB=M
source("https://www.ime.unicamp.br/~cnaber/AnaDesvTestCBM.r")



m.dados <- scan("https://www.ime.unicamp.br/~cnaber/prefauto.dat",what=list(0,0,0,0))
v.pref <- cbind(m.dados[[1]]) # preferência
v.ida <- cbind(m.dados[[2]]) # idade
v.gen <- cbind(m.dados[[3]])  # gênero
v.eciv <- cbind(m.dados[[4]]) # estado civil
n <- length(v.pref)

# Análise descritiva
genfac <- factor(as.factor(v.gen),levels=c("0","1"),labels=c("masculino","feminino")) # C(as.factor(tmater))
ecivfac <- factor(as.factor(v.eciv),levels=c("0","1"),labels=c("casado","solteiro")) # C(as.factor(tmater))
preffac <- factor(as.factor(v.pref),levels=c("0","1"),labels=c("japonês","americano"))

# por gênero
dprefgen <- data.frame(v.pref,genfac)
rpref1 <- ddply(dprefgen,.(genfac),summarise,media=mean(v.pref),dp=sqrt(var(v.pref)),vari=var(v.pref),cv=100*((sqrt(var(v.pref))/mean(v.pref))),n=length(v.pref))

# por estado civil
dprefeciv <- data.frame(v.pref,ecivfac)
rpref2 <- ddply(dprefeciv,.(ecivfac),summarise,media=mean(v.pref),dp=sqrt(var(v.pref)),vari=var(v.pref),cv=100*((sqrt(var(v.pref))/mean(v.pref))),n=length(v.pref))
xtable(rpref1)
xtable(rpref2)

# por gênero x estado civil
dpref <- data.frame(v.pref,genfac,ecivfac)
rpref3 <- ddply(dpref,.(genfac,ecivfac),summarise,media=mean(v.pref),dp=sqrt(var(v.pref)),vari=var(v.pref),cv=100*((sqrt(var(v.pref))/mean(v.pref))),n=length(v.pref))
xtable(rpref3)

#options(digits=2)
#options(op)     # reset (all) initial options
#options("digits") 


# tabelas de contingência
t1 <- table(preffac,genfac)
t2 <- table(preffac,ecivfac)
ta1 <- 100*(t(t1)/(apply(t1,2,sum)))
ta2 <- 100*(t(t2)/(apply(t2,2,sum)))
xtable(rbind(ta1,ta2))
#

# Gráficos das proporções
plot(c(ta1[1,],ta1[2,]),axes=FALSE,ylab="proporções",xlab="preferência por gênero",cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19)
axis(1,c(1,2,3,4),labels=c("masculino-japonês","masculino-americano","feminino-japonês","feminino-americano"))
#axis(2,at=c(round(min(ta1),2),round(max(ta1),2)))
axis(2,at=seq(round(min(ta1),2),round(max(ta1),2),2))

plot(c(ta2[1,],ta2[2,]),axes=FALSE,ylab="proporções",xlab="preferência por estado civil",cex=1.2,cex.axis=1.2,cex.lab=1.2,pch=19)
axis(1,c(1,2,3,4),labels=c("casado-japonês","casado-americano","solteiro-japonês","solteiro-americano"))
#axis(2,at=c(round(min(ta1),2),round(max(ta1),2)))
axis(2,at=seq(round(min(ta2),2),round(max(ta2),2),2))


# Gráfico de perfis
medias <- rpref3$media
mediasa <- medias
dp <- rpref3$dp
vn <- rpref3$n

ez <- qnorm(0.975)
plot(medias[1:2],axes=FALSE,ylim=c(0.15,0.66),cex.lab=1.5,xlab="estado civil",ylab="proporção de preferência por carros americanos")
axis(2,cex.axis=1.2)
axis(1,1:2,c("casado","solteiro"),cex.axis=1.2)
plotCI(medias[1:2],liw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[1:2]),uiw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[1:2]),pch=19,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=4,cex=1.2)
lines(medias[1:2],lwd=2,col=4)
plotCI(medias[3:4],liw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[3:4]),uiw=ez*sqrt(medias[1:2]*(1-medias[1:2]))/sqrt(vn[3:4]),pch=23,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=2,pt.bg=2,cex=1.2)
lines(medias[3:4],col=2,lwd=2)
legend(1,0.35,col=c(4,2),lwd=c(2,2),pch=c(19,23),pt.bg=c(2,2),legend=c("masculino","feminino"),bty="n",cex=1.5)



#################################
#################################
# Análise inferencial sem a idade


# Ajuste dos modelos completos (comparação dos modelos)

## logito
fit.model <- result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac,family=binomial(link="logit"))
summary(fit.model)
### OBS > genfacfeminino  = 0.1854, significa que ao mudar do sexo M pro F, vc tem um aumento de 0.1854 na probabilidade estimada de preferencia do carro japones em relação ao carro americano
###     > nota-se que nenhuma das variaveis são significativas, o que pode estar indicando que o modelo esta superparametrizado
###     > genfacfeminino:ecivfacsolteiro = representa o cruzamento da diferença quando eu tenho mudanças dentro de cada grupo; se a mudança entre casados e solteiros for a mesma,
###       independente do sexo, significa que não tem necessidade dessa iteração, ou seja, eu to assumindo que a mudança de preferência, por estado civil, não muda pros diferentes sexos,
###       ou seja, pra M e pra F, ela se mantem equivalente; a iteração é pra mostrar que essa relação muda entre os sexos, ou seja, a queda pela preferencia por carro japones, do grupo 
###       de casado em relação ao grupo de solteiro, é diferente no sexo M em relação ao sexo F

diagBern(fit.model)
envelBern(fit.model,"logit")

rebeta <- (summary(result))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4]) # incluindo os IC's
xtable(rebetaLL)

confint(fit.model) # outra maneira, mais direta, de calcular os IC's

AICLL <- AIC(fit.model)
BICLL <- BIC(fit.model)

mrebeta <- rbind(rebetaLL)
pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit

## probito
fit.model <- result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac,family=binomial(link="probit"))
summary(fit.model)

diagBern(fit.model)
envelBern(fit.model,"probit")

rebeta <- (summary(result))$coef
rebetaLP <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLP)

AICLP <- AIC(fit.model)
BICLP <- BIC(fit.model)

mrebeta <- rbind(rebetaLP)
predLP <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLP <- predLP$fit

## cauchito
fit.model <- result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac,family=binomial(link="cauchit"))
summary(fit.model)

diagBern(fit.model)
envelBern(fit.model,"cauchit")

rebeta <- (summary(result))$coef
rebetaLC <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLC)

AICLC <- AIC(fit.model)
BICLC <- BIC(fit.model)

mrebeta <- rbind(rebetaLC)
predLC <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLC <- predLC$fit

## cloglog
fit.model <- result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac,family=binomial(link="cloglog"))
summary(fit.model)

diagBern(fit.model)
envelBern(fit.model,"cloglog")

rebeta <- (summary(result))$coef
rebetaLCLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLCLL)

AICLCLL <-AIC(fit.model)
BICLCLL <- BIC(fit.model)

mrebeta <- rbind(rebetaLCLL)
predLCLL <-predict(fit.model,type=c("response"),se.fit = TRUE)
mupredLCLL <- predLCLL$fit


# AIC e BIC
AICBIC <- rbind(cbind(AICLL,AICLP,AICLC,AICLL),cbind(BICLL,BICLP,BICLC,BICLL))
colnames(AICBIC) <- c("logito","probito","cauchito","cloglog")
rownames(AICBIC) <- c("AIC","BIC")
xtable(t(AICBIC))

# Desvio absolutos médios
damLL <- mean(abs(v.pref-mupred))
damLP <- mean(abs(v.pref-mupredLP))
damLC <- mean(abs(v.pref-mupredLC))
damLCLL <- mean(abs(v.pref-mupredLCLL))

AICBICDAM <- rbind(cbind(AICLL,AICLP,AICLC,AICLCLL),cbind(BICLL,BICLP,BICLC,BICLCLL),cbind(damLL,damLP,damLC,damLCLL))
colnames(AICBICDAM) <- c("logito","probito","cauchito","cloglog")
rownames(AICBICDAM) <- c("AIC","BIC","DAM")
xtable(t(AICBICDAM))



####################
# modelo logito

#sem interação
fit.model <- result1 <- glm(v.pref~genfac+ecivfac,family=binomial(link="logit"))
summary(fit.model)
### OBS > componente de estado civil passa a ser significativa

diagBern(fit.model)
envelBern(fit.model,"logit")

rebeta <- (summary(result1))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)

mrebeta <- rbind(rebetaLL)
pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit


#somente com estado civil
fit.model <- result2 <- glm(v.pref~ecivfac,family=binomial(link="logit"))
summary(fit.model)

diagBern(fit.model)
envelBern(fit.model,"logit")

rebeta <- (summary(result2))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)

mrebeta <- rbind(rebetaLL)
pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit


# Perfis ajustados
rpmedias <- predict(fit.model,type=c("response"),se.fit = TRUE)
medias <- unique(rpmedias$fit)
epmedias <- unique(rpmedias$se.fit)


plot(medias[1:2],axes=FALSE,ylim=c(0.25,0.60),cex.lab=1.5,xlab="estado civil",ylab="proporção de preferência por carros americanos")
axis(2,cex.axis=1.2)
axis(1,1:2,c("casado","solteiro"),cex.axis=1.2)
plotCI(medias[1:2],liw=ez*epmedias[1:2],uiw=ez*epmedias[1:2],pch=19,add=TRUE,cex.lab=1.2,slty=1,lwd=2,col=1,cex=1.2)
#plotCI(medias[1:2],liw=1.96*medias[1:2]*dp[1:2]/sqrt(vn[1:2]),uiw=1.96*medias[1:2]*dp[1:2]/sqrt(vn[1:2]),pch=19,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=4,cex=1.2)
lines(medias[1:2],lwd=2,col=1)
lines(mediasa[1:2],type="p",pch=15,col=2,cex=1.2)
lines(mediasa[3:4],type="p",pch=16,col=3,cex=1.2)
#plotCI(medias[3:4],liw=ez*epmedias[3:4],uiw=ez*epmedias[3:4],pch=23,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=2,pt.bg=2,cex=1.2)
#plotCI(medias[3:4],liw=1.96*medias[3:4]*dp[3:4]/sqrt(vn[3:4]),uiw=1.96*medias[3:4]*dp[3:4]/sqrt(vn[3:4]),pch=23,add=TRUE,cex.lab=1.5,slty=1,lwd=2,col=2,pt.bg=2,cex=1.2)
#lines(medias[3:4],col=2,lwd=2)
#legend(1.5,0.50,col=c(4,2),lwd=c(2,2),pch=c(19,23),pt.bg=c(2,2),legend=c("masculino","feminino"),bty="n",cex=1.5)
#legend(1.5,0.50,col=c(1),lwd=c(2),pch=c(19),pt.bg=c(2),legend=c("masculino/feminino"),bty="n",cex=1.5)
legend(1,0.35,col=c(1,2,3),lwd=c(2,0,0),lty=c(1,0,0),pch=c(19,15,16),pt.bg=c(2,2,2),legend=c("predito - masculino/feminino","observado - masculino","observado - feminino"),bty="n",cex=1.4)


mmedias <- cbind(c(medias,medias),c(epmedias,epmedias),c(medias,medias)-ez*c(epmedias,epmedias),c(medias,medias)+ez*c(epmedias,epmedias))
mmedias <- round(100*mmedias,2)
alabel1 <- c("Casado","Solteiro","Casado","Solteiro")
alabel2 <- c("Masculino","Masculino","Feminino","Feminino")
xtable(cbind(alabel1,alabel2,mmedias))


#modelo completo
fit.model<-result<- glm(v.pref~genfac+ecivfac+genfac*ecivfac,family=binomial(link="logit"))

#modelo somente com o intercepto
fit.model0<-result0<- glm(v.pref~1,family=binomial(link="logit"))

# Seleção automatizada

step(fit.model0,scope=list(lower=result0,upper=result),direction=c("forward"))
step(fit.model,direction=c("backward"))
step(fit.model0,scope=list(upper=result),direction=c("both"))



################################
################################
################################
# Análise contemplando a idade


# boxplots
par(mfrow=c(2,2))
boxplot(v.ida~v.pref,cex=1.2,cex.lab=1.2,cex.axis=1.2,xlab="preferência", ylab="idade",names=c("japonês","americano"))
boxplot(v.ida~v.pref*genfac,cex=1.2,cex.lab=1.2,cex.axis=1.2,xlab="preferência", ylab="idade",names=c("M - japonês","M - americano","F - japonês","F - americano"))
boxplot(v.ida~v.pref*ecivfac,cex=1.2,cex.lab=1.2,cex.axis=1.2,xlab="preferência", ylab="idade",names=c("C - japonês","C - americano","S - japonês","S - americano"))
#boxplot(v.ida~v.pref*genfac*ecivfac,cex=1.2,cex.lab=1.2,cex.axis=1.2,xlab="preferência", ylab="idade",names=c("MC/jap.","MC/amer.","FC/jap.","FC/amer.","MS/jap.","MS/amer.","FS/jap.","FS/amer."))
boxplot(v.ida~v.pref*genfac*ecivfac,cex=1.2,cex.lab=1.2,cex.axis=1.2,xlab="preferência", ylab="idade",names=c("MC-J","MC-A","FC-J","FC-A","MS-J","MS-A","FS-J","FS-A"))


# Gráficos de dispersão
aux <- data.frame(cbind(v.pref,as.factor(v.ida),genfac))
obsfc <- c(v.pref[genfac=="feminino" & ecivfac=="casado"])
obsfs <- c(v.pref[genfac=="feminino" & ecivfac=="solteiro"])
obsmc <- c(v.pref[genfac=="masculino" & ecivfac=="casado"])
obsms <- c(v.pref[genfac=="masculino" & ecivfac=="solteiro"])

par(mfrow=c(2,2))
plot(v.ida[genfac=="feminino" & ecivfac=="casado"],obsfc,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="feminino-casado")
plot(v.ida[genfac=="feminino" & ecivfac=="solteiro"],obsfs,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="feminino-solteiro")
plot(v.ida[genfac=="masculino" & ecivfac=="casado"],obsmc,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="masculino-casado")
plot(v.ida[genfac=="masculino" & ecivfac=="solteiro"],obsms,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="maculino-solteiro")


#propofc<-c(by(aux[genfac=="feminino" & ecivfac=="casado",1],aux[genfac=="feminino"& ecivfac=="casado",2],mean))
#propofs<-c(by(aux[genfac=="feminino" & ecivfac=="solteiro",1],aux[genfac=="feminino"& ecivfac=="solteiro",2],mean))
#propomc<-c(by(aux[genfac=="masculino"& ecivfac=="casado",1],aux[genfac=="masculino"& ecivfac=="casado",2],mean))
#propoms<-c(by(aux[genfac=="masculino"& ecivfac=="solteiro",1],aux[genfac=="masculino"& ecivfac=="solteiro",2],mean))

#par(mfrow=c(2,2))
#plot(aux[genfac=="feminino"& ecivfac=="casado",2],propofc,xlab="idade",ylab="proporção de preferência por automóveis americanos",main="feminino-casado")


par(mfrow=c(2,2))
cfc <- table(v.ida[genfac=="feminino" & ecivfac=="casado"])
cfc1 <- table(v.pref[genfac=="feminino" & ecivfac=="casado"],v.ida[genfac=="feminino" & ecivfac=="casado"])
plot(as.vector(count(cfc)[,1]),as.numeric(cfc1[2,]/cfc),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="feminino-casado")

cfs <- table(v.ida[genfac=="feminino" & ecivfac=="solteiro"])
cfs1 <- table(v.pref[genfac=="feminino" & ecivfac=="solteiro"],v.ida[genfac=="feminino" & ecivfac=="solteiro"])
plot(as.vector(count(cfs)[,1]),as.numeric(cfs1[2,]/cfs),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="feminino-solteiro")

cmc <- table(v.ida[genfac=="masculino" & ecivfac=="casado"])
cmc1 <- table(v.pref[genfac=="masculino" & ecivfac=="casado"],v.ida[genfac=="masculino" & ecivfac=="casado"])
plot(as.vector(count(cmc)[,1]),as.numeric(cmc1[2,]/cmc),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="masculino-casado")

cms<-table(v.ida[genfac=="masculino" & ecivfac=="solteiro"])
cms1 <- table(v.pref[genfac=="masculino" & ecivfac=="solteiro"],v.ida[genfac=="masculino" & ecivfac=="solteiro"])
plot(as.vector(count(cms)[,1]),as.numeric(cms1[2,]/cms),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="masculino-solteiro")



#####################
#####################
# Análise inferencial


# Modelo com a idade com um coeficiente para cada grupo.
mida <- mean(v.ida) # idade media
auxida <- v.ida-mida # variavel corrigida da idade media
auxida2 <- auxida^2

#fit.model <- result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac+auxida,family=binomial(link="logit"))
fit.model <-result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac+auxida*genfac*ecivfac,family=binomial(link="logit"))
summary(fit.model)

diagBern(fit.model)
envelBern(fit.model,"logit")

rebeta <- (summary(result))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)

AIC(fit.model)
BIC(fit.model)

mrebeta <- rbind(rebetaLL)


# Único coeficiente para a idade ; sem fazer iteração entre idade e os demais componentes
fit.model0 <- result <- glm(v.pref~genfac+ecivfac+genfac*ecivfac+auxida,family=binomial(link="logit"))
summary(fit.model0)
### OBS > a idade entra como fator positivo, ou seja, para cada incremento da idade em relação a idade media, há uma tendência de aumentar em 0.05109
###       a probabilidade de preferência  pelo carro americano em relação ao carro japones

diagBern(fit.model0)
envelBern(fit.model0,"logit")

rebeta <- (summary(result))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)

AIC(fit.model0)
BIC(fit.model0)

mrebeta <- rbind(rebetaLL)


# Análise do desvio
anadesv(fit.model0,fit.model) # não ha evidencias para rejeitar h0, entao optamos por considerar esses coeficientes de iteração iguais a zero, entao optamos pelo modelo mais parcimonioso


# Único coeficiente para a idade sem interação gênero x estado civil
fit.model <- result <- glm(v.pref~genfac+ecivfac+auxida,family=binomial(link="logit"))
summary(fit.model)

diagBern(fit.model)
envelBern(fit.model,"logit")

rebeta <- (summary(result))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)

AICLL <- AIC(fit.model)
BICLL <- BIC(fit.model)

mrebeta <- rbind(rebetaLL)
pred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit


# Único coeficiente para a idade sem interação gênero x estado civil e sem o fator gênero
fit.model <- result <- glm(v.pref~ecivfac+auxida,family=binomial(link="logit"))
summary(fit.model)
### OBS > ecivfacsolteiro = -0.52578: quando passa de casado pra solteiro, tem-se uma queda na preferencia por carros americanos, ou seja, o grupo de casados tem 
###       uma probabilidade maior da preferencia de carro americano em relação ao grupo de solteiros
###     > auxida = 0.04959: a medida que a idade aumenta em relação a sua media, tem-se uma aumento na probabilidade de preferencia pelo carro americano em relação ao carro japones

diagBern(fit.model)
envelBern(fit.model,"logit")

rebeta <- (summary(result))$coef
rebetaLL <- cbind(rebeta[,1],rebeta[,2],rebeta[,1]-ez*rebeta[,2],rebeta[,1]+ez*rebeta[,2],rebeta[,3],rebeta[,4])
xtable(rebetaLL)

AICLL <- AIC(fit.model)
BICLL <- BIC(fit.model)

mrebeta <- rbind(rebetaLL)


# Médias preditas pelo modelo reduzido
fit.model <- result <- glm(v.pref~ecivfac+auxida,family=binomial(link="logit"))

rmupred <- predict(fit.model,type=c("response"),se.fit = TRUE)
mupred <- pred$fit
semupred <- pred$se.fit
LIIC <- mupred-ez*semupred
LSIC <- mupred+ez*semupred


par(mfrow=c(2,2))
cfc <- table(v.ida[genfac=="feminino" & ecivfac=="casado"])
cfc1 <- table(v.pref[genfac=="feminino" & ecivfac=="casado"],v.ida[genfac=="feminino" & ecivfac=="casado"])
cfc2 <- as.numeric(rownames(as.matrix(table(mupred[genfac=="feminino" & ecivfac=="casado"],v.ida[genfac=="feminino" & ecivfac=="casado"]))))
cfc3 <- as.numeric(rownames(as.matrix(table(LIIC[genfac=="feminino" & ecivfac=="casado"],v.ida[genfac=="feminino" & ecivfac=="casado"]))))
cfc4 <- as.numeric(rownames(as.matrix(table(LSIC[genfac=="feminino" & ecivfac=="casado"],v.ida[genfac=="feminino" & ecivfac=="casado"]))))
plot(as.vector(count(cfc)[,1]),as.numeric(cfc1[2,]/cfc),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="feminino-casado",ylim=c(min(LIIC),max(LSIC)))
text(as.vector(count(cfc)[,1]),as.numeric(cfc1[2,]/cfc-0.05),as.vector(count(cfc)[,2]),cex=1.2)
plotCI(as.numeric(as.vector(count(cfc)[,1])),cfc2,li=cfc3,ui=cfc4,add=TRUE,col=2,pch=19,cex=1.2)


cfs <- table(v.ida[genfac=="feminino" & ecivfac=="solteiro"])
cfs1 <- table(v.pref[genfac=="feminino" & ecivfac=="solteiro"],v.ida[genfac=="feminino" & ecivfac=="solteiro"])
cfs2 <- as.numeric(rownames(as.matrix(table(mupred[genfac=="feminino" & ecivfac=="solteiro"],v.ida[genfac=="feminino" & ecivfac=="solteiro"]))))
cfs3 <- as.numeric(rownames(as.matrix(table(LIIC[genfac=="feminino" & ecivfac=="solteiro"],v.ida[genfac=="feminino" & ecivfac=="solteiro"]))))
cfs4 <- as.numeric(rownames(as.matrix(table(LSIC[genfac=="feminino" & ecivfac=="solteiro"],v.ida[genfac=="feminino" & ecivfac=="solteiro"]))))
plot(as.vector(count(cfs)[,1]),as.numeric(cfs1[2,]/cfs),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="feminino-solteiro",ylim=c(min(LIIC)-0.05,max(LSIC)))
text(as.vector(count(cfs)[,1]),as.numeric(cfs1[2,]/cfs-0.05),as.vector(count(cfs)[,2]),cex=1.2)
plotCI(as.numeric(as.vector(count(cfs)[,1])),cfs2,li=cfs3,ui=cfs4,add=TRUE,col=2,pch=19,cex=1.2)


cmc <- table(v.ida[genfac=="masculino" & ecivfac=="casado"])
cmc1 <- table(v.pref[genfac=="masculino" & ecivfac=="casado"],v.ida[genfac=="masculino" & ecivfac=="casado"])
cmc2 <- as.numeric(rownames(as.matrix(table(mupred[genfac=="masculino" & ecivfac=="casado"],v.ida[genfac=="masculino" & ecivfac=="casado"]))))
cmc3 <- as.numeric(rownames(as.matrix(table(LIIC[genfac=="masculino" & ecivfac=="casado"],v.ida[genfac=="masculino" & ecivfac=="casado"]))))
cmc4 <- as.numeric(rownames(as.matrix(table(LSIC[genfac=="masculino" & ecivfac=="casado"],v.ida[genfac=="masculino" & ecivfac=="casado"]))))
plot(as.vector(count(cmc)[,1]),as.numeric(cmc1[2,]/cmc),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="masculino-casado",ylim=c(min(LIIC)-0.05,max(LSIC)))
text(as.vector(count(cmc)[,1]),as.numeric(cmc1[2,]/cmc-0.05),as.vector(count(cmc)[,2]),cex=1.2)
plotCI(as.numeric(as.vector(count(cmc)[,1])),cmc2,li=cmc3,ui=cmc4,add=TRUE,col=2,pch=19,cex=1.2)


cms <- table(v.ida[genfac=="masculino" & ecivfac=="solteiro"])
cms1 <- table(v.pref[genfac=="masculino" & ecivfac=="solteiro"],v.ida[genfac=="masculino" & ecivfac=="solteiro"])
cms2 <- as.numeric(rownames(as.matrix(table(mupred[genfac=="masculino" & ecivfac=="solteiro"],v.ida[genfac=="masculino" & ecivfac=="solteiro"]))))
cms3 <- as.numeric(rownames(as.matrix(table(LIIC[genfac=="masculino" & ecivfac=="solteiro"],v.ida[genfac=="masculino" & ecivfac=="solteiro"]))))
cms4 <- as.numeric(rownames(as.matrix(table(LSIC[genfac=="masculino" & ecivfac=="solteiro"],v.ida[genfac=="masculino" & ecivfac=="solteiro"]))))
plot(as.vector(count(cms)[,1]),as.numeric(cms1[2,]/cms),pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,xlab="idade",ylab="preferência por automóveis americanos",main="masculino-solteiro",ylim=c(min(LIIC)-0.05,max(LSIC)))
text(as.vector(count(cms)[,1]),as.numeric(cms1[2,]/cms-0.05,as.vector(count(cms)[,2])),cex=1.2)
plotCI(as.numeric(as.vector(count(cms)[,1])),cms2,li=cms3,ui=cms4,add=TRUE,col=2,pch=19,cex=1.2)


# Médias preditas
#table(v.ida[genfac=="feminino"])
#table(v.ida[genfac=="masculino"])

#mobsF <- c(mean(v.pref[genfac=="feminino" & v.ida==23]),mean(v.pref[genfac=="feminino" & v.ida==24]),
#mean(v.pref[genfac=="feminino" & v.ida==25]),mean(v.pref[genfac=="feminino" & v.ida==26]),
#mean(v.pref[genfac=="feminino" & v.ida==27]),mean(v.pref[genfac=="feminino" & v.ida==28]),
#mean(v.pref[genfac=="feminino" & v.ida==29]),mean(v.pref[genfac=="feminino" & v.ida==30]),
#mean(v.pref[genfac=="feminino" & v.ida==31]),mean(v.pref[genfac=="feminino" & v.ida==32]),
#mean(v.pref[genfac=="feminino" & v.ida==33]))

#mpredF <- c(mean(vYpred[genfac=="feminino" & v.ida==23]),mean(vYpred[genfac=="feminino" & v.ida==24]),
#mean(vYpred[genfac=="feminino" & v.ida==25]),mean(vYpred[genfac=="feminino" & v.ida==26]),
#mean(vYpred[genfac=="feminino" & v.ida==27]),mean(vYpred[genfac=="feminino" & v.ida==28]),
#mean(vYpred[genfac=="feminino" & v.ida==29]),mean(vYpred[genfac=="feminino" & v.ida==30]),
#mean(vYpred[genfac=="feminino" & v.ida==31]),mean(vYpred[genfac=="feminino" & v.ida==32]),
#mean(vYpred[genfac=="feminino" & v.ida==33]))
#plot(mobsF)
#lines(mpredF,type="p",pch=19)


#mobsm <- c(mean(v.pref[genfac=="feminino" & v.ida==23]),mean(v.pref[genfac=="feminino" & v.ida==24]),
#mean(v.pref[genfac=="feminino" & v.ida==25]),mean(v.pref[genfac=="feminino" & v.ida==26]),
#mean(v.pref[genfac=="feminino" & v.ida==27]),mean(v.pref[genfac=="feminino" & v.ida==28]),
#mean(v.pref[genfac=="feminino" & v.ida==29]),mean(v.pref[genfac=="feminino" & v.ida==30]),
#mean(v.pref[genfac=="feminino" & v.ida==31]),mean(v.pref[genfac=="feminino" & v.ida==32]),
#mean(v.pref[genfac=="feminino" & v.ida==33]),mean(v.pref[genfac=="feminino" & v.ida==34]),
#mean(v.pref[genfac=="feminino" & v.ida==35]),mean(v.pref[genfac=="feminino" & v.ida==38]))
#plot(mean(v.pref[genfac=="feminino" & v.ida==23]))
#lines(mean(v.pref[genfac=="feminino" & v.ida==23]))