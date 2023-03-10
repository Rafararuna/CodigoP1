########################################
########################################
# Dados de pagamento de seguros
########################################
########################################

# Pacotes
library(xtable)
library(TeachingDemos)
library(plotrix)
library(plyr)
library(e1071)
library(MASS) 


# Diagn?stico e envelope para o modelo gama
source("https://www.ime.unicamp.br/~cnaber/diag_gama.r")
source("https://www.ime.unicamp.br/~cnaber/envel_gama.r")
# Fun??es auxiliares para o modelo gama
source("https://www.ime.unicamp.br/~cnaber/AuxModGama.r")
# Fun??esde An?lise do Desvio e do teste CB=M
source("https://www.ime.unicamp.br/~cnaber/AnaDesvTestCBM.r")


# Carregando os dados
m.dados <- scan("https://www.ime.unicamp.br/~cnaber/insurance.prn",what=list(0,0,0,0))


# Gerando e/ou nomeando as vari?veis
vpago <- cbind(m.dados[[1]]) 
optime <- cbind(m.dados[[4]])
legrep <- cbind(m.dados[[2]])
legrepfac <- factor(legrep, c("0","1"))
legrepfac <- revalue(legrepfac, c("0" = "nao", "1" = "sim"))
n <- length(vpago)


# Imprimindo os dados em forma de tabela
mdados <- cbind(vpago,optime,legrep)
xtable(mdados)



####################
####################
# An?lise descritiva
####################
####################


par(mfrow=c(1,2))
plot(optime[legrepfac=="nao"],vpago[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,
     xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",
     main="Sem representa??o legal",pch=19) # no grafico ? esquerda tamb?m notamos uma tend?ncia de acrescimo 
                                            # nota-se uma alta dispers?o 
                                            # no final parece que os pontos caem (os maiores valores est?o na regi?o central), ent?o tem-se um tend?ncia de crescimento e dps queda
plot(optime[legrepfac=="sim"],vpago[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,
     xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",
     main="Com representa??o legal",pch=19) # nota-se uma concentra??o dos dados na parte inferior do gr?fico, bem como uma tend?ncia de acr?scimo, ou seja, a medida que o tempo operacional aumenta tem uma tend?ncia de o valor pago ser maior
                                            # nota-se pontos com valores muito altos

# boxplot (rela??o de cada variavel numerica com a variavel categorica)
par(mfrow=c(1,2))
boxplot(vpago~legrepfac,names=c("n?o","sim"),xlab="representa??o legal",ylab="valor pago ao segurado",
        cex=1.2,cex.axis=1.2,cex.lab=1.2) # nota-se um aumento no grupo "sim"
                                          # em rela??o a dispers?o, parece ser um pouco menor no grupo "sim"
boxplot(optime~legrepfac,names=c("n?o","sim"),xlab="representa??o legal",
        ylab="tempo operacional para pagamento do seguro",cex=1.2,cex.axis=1.2,cex.lab=1.2) # n?o ha muita diferena em rela??o em termos da distribui??o n?o h? mtas mudan?as, em termos de valores m?nimos e maximos
                                                                                            # nota-se um deslocamento das medimas resumo, terceiro e primeiro quantil, mediana, essas medidas estao acima no grupo "sim", porem mto sutil essa diferen?a

table(legrepfac)
dataseg <- data.frame(vpago,optime,legrepfac)

# Medidas resumo
cseg1 <- ddply(dataseg,.(legrepfac),summarise,media=mean(vpago),dp=sqrt(var(vpago)),
               vari=var(vpago),cv=100*((sqrt(var(vpago))/mean(vpago))),
               ca=skewness(vpago),minimo=min(vpago),maximo=max(vpago)) # medidas resumo da variavel vpago por grupo
cseg2 <- ddply(dataseg,.(legrepfac),summarise,media=mean(optime),dp=sqrt(var(optime)),
               vari=var(optime),cv=100*((sqrt(var(optime))/mean(optime))),
               ca=skewness(optime),minimo=min(optime),maximo=max(optime)) # medidas resumo da variavel optime por grupo

xtable(rbind(cseg1,cseg2))
xtable(cbind(t(cseg1),t(cseg2)))



#####################
#####################
# An?lise Inferencial
#####################
#####################


########################
########################
# Modelo Gama (link log)

optime2 <- optime^2 # valor da variavel optime ao quadrado
optimec <- optime-0.1 # valor da variavel optime corrigido
optimec2 <- (optime-0.1)^2 # valor da variavel optime ao quadrado corrigido 
#optime2 <- ifelse(legrepfac=="nao",optime2,0)
fit.model <- glm(vpago~legrepfac+optimec+optimec2+legrepfac*optimec+legrepfac*optimec2,family=Gamma(link="log"))
#fit.model <- glm(vpago~legrepfac+optime+optime2+legrepfac*optime+legrepfac*optime2,family=Gamma(link="log")) # modelo sem corrigir o optime, ou seja, sem subtrair 0.1
#fit.model <- glm(vpago~legrepfac+optime+optime2+legrepfac*optime,family=Gamma(link="log")) # modelo sem a intera??o de segunda ordem e com os optimes sem corre??o
summary(fit.model)


# Usando a estimativa de MV de phi
summary(fit.model,dispersion = gamma.dispersion(fit.model))
xtable(summary(fit.model,dispersion = gamma.dispersion(fit.model)))
p <-ncol(model.matrix(fit.model)) 


# Estimativa do par?metro phi
# Default R
ephi <- 1/summary(fit.model)$dispersion

# M?todo dos momentos
ro <- resid(fit.model,type="response")
ephiM <-(n-p)/sum((ro/(fitted(fit.model)))^ 2)


# MV
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 
sephiMV <- resultphiMV$SE

ez <- qnorm(0.975)
ICphi <- cbind(ephiMV-ez*sephiMV,ephiMV+ez*sephiMV)


# N?o consistente
ephiNC <- (n-p)/deviance(fit.model)


# Estimativas dos par?metros
resultmg <- coef(summary(fit.model,dispersion = gamma.dispersion(fit.model)))
ez <- qnorm(0.975)
ICbeta <-cbind(resultmg[,1]-ez*resultmg[,2],resultmg[,1]+ez*resultmg[,2])
resultmg <- cbind(resultmg[,1],resultmg[,2],ICbeta,resultmg[,3],resultmg[,4])
resultmg <- rbind(resultmg,cbind(c(ephiMV),c(sephiMV),ICphi,0,0))
xtable(resultmg)


# AIC e BIC
AICBICgama(vpago,fitted(fit.model),ephiMV,n,p+1)


# Deviance do modelo
desv <- deviance(fit.model)*ephiMV # desvio escalonado ; se fosse pra calcula o n?o escalonado era s? n?o multiplicar por ephiMV
pvalordesv <- 1-pchisq(desv,df=n-p) # p-valor baixo, rejeita H0, ou seja, o modelo n?o est? adequado, est? rejeitando o ajuste do modelo


# Gr?ficos de diagn?stico
diaggama(fit.model,1) # em ambos os gr?ficos de cima nota-se a presen?a de muitos outliers
                      # no grafico dos residuso pelo valor ajustado nota-se uma mudan?a na variabilidade (heterocedasticidade), que come?a co =m alta dispers?o e vai diminuindo
envelgama(fit.model,"log",1) # nota-se que nesse gr?fico tbm h? um problema no ajuste, visto que muitos pontos est?o fora da banda


# M?dias preditas pelo modelo
resmupred <- predict(fit.model,type = c("response"),se.fit=TRUE)
mupred1 <- (resmupred$fit)
epmupred1 <- (resmupred$se.fit)

par(mfrow=c(1,1))
plot(optime[legrepfac=="nao"],mupred1[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Valor predito por tempo operacional",pch=19,type="l",lty=1,ylim=c(min(mupred1),max(mupred1)),lwd=2)
lines(optime[legrepfac=="sim"],mupred1[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",pch=19,lty=2,lwd=2)
legend(1,13000,c("sem representa??o legal","com representa??o legal"),col=c(1,1),lty=c(1,2),lwd=c(3,3),bty="n",cex=1.2)
# OBS > notamos que o grupo "com representa??o legal" aparenta ter um comportamento linear
#     > notamos que o grupo "sem representa??o legal" aparenta ter um comportamento quadr?tico

par(mfrow=c(1,2))
plot(optime[legrepfac=="nao"],vpago[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Sem representa??o legal",pch=19,ylim=c(-1000,40000))
#plotCI(optime[legrepfac=="nao"],mupred1[legrepfac=="nao"],uiw=1.96*mupred1[legrepfac=="nao"],liw=1.96*mupred1[legrepfac=="nao"],add=TRUE,col=2)
lines(optime[legrepfac=="nao"],mupred1[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Sem representa??o legal",pch=19,type="l",lwd=3,col=2)
plot(optime[legrepfac=="sim"],vpago[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Com representa??o legal",pch=19)
lines(optime[legrepfac=="sim"],mupred1[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Com representa??o legal",pch=19,type="l",lwd=3,col=2)


# O objetivo ? obter os IC assint?ticos
par(mfrow=c(2,2))
plot(optime[legrepfac=="nao"],vpago[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Sem representa??o legal (IC assint?tico)",pch=19,ylim=c(-10000,40000))
plotCI(optime[legrepfac=="nao"],mupred1[legrepfac=="nao"],uiw=1.96*epmupred1[legrepfac=="nao"],liw=1.96*epmupred1[legrepfac=="nao"],add=TRUE,col=2,pch=19)
#lines(optime[legrepfac=="nao"],mupred1[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Sem representa??olegal",pch=19,type="l",lwd=3,col=2)
plot(optime[legrepfac=="sim"],vpago[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Com representa??o legal (IC assint?tico)",pch=19,ylim=c(-20000,120000))
plotCI(optime[legrepfac=="sim"],mupred1[legrepfac=="sim"],uiw=1.96*epmupred1[legrepfac=="sim"],liw=1.96*epmupred1[legrepfac=="sim"],add=TRUE,col=2,pch=19)
#lines(optime[legrepfac=="sim"],mupred1[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Com representa??olegal",pch=19,type="l",lwd=3,col=2)


# Outra forma de calcular os IC assint?ticos ? por meio de simula??o
# Predi??o via simula??o (o obketivo ? obter um Intervalor de confian?a empirico)
X <- model.matrix(fit.model)
resp <- rgamma(n,ephiMV)
resp <- (fitted(fit.model)/ephiMV)*resp
result <- glm(resp ~ X,family = Gamma(link = "log"))

mupred <- as.vector(predict(result,type = c("response"),se.fit=TRUE)$fit)

for (k in 1:99)
{
  resp <- rgamma(n,ephiMV)
  resp <- (fitted(fit.model)/ephiMV)*resp
  result<-glm(resp ~ X,family=Gamma(link="log"))
  mupred <- rbind(mupred,as.vector(predict(result,type = c("response"),se.fit=TRUE)$fit))
  #mupred <-as.matrix(mupred)
}

mmupred <- apply(mupred,2,mean) # media predita
liicmupred <-apply(mupred,2,quantile,0.025) # limite inferior predito
lsicmupred <-apply(mupred,2,quantile,0.975) # limite superior predito

par(mfrow=c(1,2))
plot(optime[legrepfac=="nao"],vpago[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Sem representa??o legal (IC via predi??o)",pch=19,ylim=c(-1000,40000))
plotCI(optime[legrepfac=="nao"],mmupred[legrepfac=="nao"],li=liicmupred[legrepfac=="nao"],ui=lsicmupred[legrepfac=="nao"],add=TRUE,col=2,pch=19)
#lines(optime[legrepfac=="nao"],mmupred[legrepfac=="nao"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Sem presenta??????o legal",pch=19,type="l",lwd=3,col=2)
plot(optime[legrepfac=="sim"],vpago[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Com representa??o legal (IC via predi??o)",pch=19)
plotCI(optime[legrepfac=="sim"],mmupred[legrepfac=="sim"],li=liicmupred[legrepfac=="sim"],ui=lsicmupred[legrepfac=="sim"],add=TRUE,col=2,pch=19)
#lines(optime[legrepfac=="sim"],mmupred[legrepfac=="sim"],cex=1.2,cex.lab=1.2,cex.main=1.2,xlab="tempo operacional para pagamento do seguro",ylab="valor pago ao segurado",main="Com presenta??????o legal",pch=19,type="l",lwd=3,col=2)

# OBS > comparando os graficos do IC assintotica com os graficos do IC via predi??o, nota-se que n?o muitas mudan?as bruscas,
#       aparentam ter comportamento semelhante no geral



############
###########
# Quest?es extras
###########

# Como definir o IC do tempo operacional que maximiza o valor pago ?
# Como definir o IC pro valor pago m?ximo ? 

# Respostas para o grupo "sem representa??o legal":

# Valor ?timo do tempo operacional para pagamento do seguro, o qual
# retorna o valor m?ximo pago do seguro, para o grupo sem representa??o legal

emax <- 0.1-resultmg[3,1]/(2*resultmg[4,1]) # (0.1 - beta/(2gama)), onde beta = optimec e gama = optimec2
mcovbeta <- (vcov(fit.model)*ephi)/(ephiMV) # matriz de variancia e covariancia em termos de phi estimado por MV ; por isso dividiu-se por "ephiMV"
auxvar <- cbind(0,0,-1/(2*resultmg[4,1]),resultmg[3,1]/(2*resultmg[4,1]^2),0,0) # gradiente ; deriva-se emax em rela??o a todos os betas
                                                                                # primeiro derivou-se em rela??o ? aplha e alpha 2, ambos deram zero
                                                                                # depois derivou-se em rela??o ? beta e gama
                                                                                # depois derivou-se em rela??o ? beta2 e gama2, ambos deram zero
varemax <- auxvar%*%mcovbeta%*%t(auxvar) # opera??o quadr?tica ; termo quadr?tico
epvaremax <- sqrt(varemax) # erro padr?o
ICemax <- c(emax-ez*epvaremax,emax+ez*epvaremax) # Intervalo de confian?a, onde ez = 1,96

# Valor m?ximo do seguro pago para o grupo sem representa??o legal

vbeta <- cbind(coef(summary(fit.model,dispersion = gamma.dispersion(fit.model)))[,1]) # pegando os coeficientes estimados dos betas
vpagomax <- exp(vbeta[1]+(emax-0.1)*vbeta[3] + (emax-0.1)^2*vbeta[4]) 
auxep <- cbind(vpagomax,0,vpagomax*(emax-0.1),vpagomax*(emax-0.1)^2,0,0) # gradiente ; deriva-se vpagomax em rela??o a todos os betas
epvpagomax <- sqrt(auxep%*%(vcov(fit.model)*ephi/ephiMV)%*%t(auxep)) # erro padr?o
ICvpagomax <- c(vpagomax-ez*epvpagomax,vpagomax+ez*epvpagomax) # Intervalo de confian?a, onde ez = 1,96