########################################
########################################
########################################
# Dados da potencia de turbina de aviões
########################################
########################################
########################################


# Código desenvolvido pelo prof. Caio Azevedo


# Carregando pacotes:
library(xtable)
library(TeachingDemos)
library(plotrix)
library(plyr)
library(e1071)
library(MASS) 


# Diagnóstico e envelope para o modelo normal
source("https://www.ime.unicamp.br/~cnaber/diag_norm.r")
source("https://www.ime.unicamp.br/~cnaber/envel_norm.r")
# Diagnóstico e envelope para o modelo gama
source("https://www.ime.unicamp.br/~cnaber/diag_gama.r")
source("https://www.ime.unicamp.br/~cnaber/envel_gama.r")
# Funções auxiliares para o modelo gama
source("https://www.ime.unicamp.br/~cnaber/AuxModGama.r")
# Funções de Análise do Desvio e do teste CB=M
source("https://www.ime.unicamp.br/~cnaber/AnaDesvTestCBM.r")


# Carregando os dados
m.dados <- scan("https://www.ime.unicamp.br/~cnaber/Turbina.prn",what=list(0,0))


## Gerando e/ou nomeando as variáveis
veloc <- tempo <- cbind(m.dados[[2]])
grupo <- cbind(m.dados[[1]])
grupo_fac <- factor(grupo,c("1","2","3","4","5"))


# Imprimindo os dados em forma de tabela
mdados <- cbind(veloc,grupo)
xtable(mdados)



####################
####################
# Análise descritiva
####################
####################

dataturb <- data.frame(veloc, grupo_fac)


# Medidas resumo
cturb <- ddply(dataturb, .(grupo_fac), summarise, media = mean(veloc), 
               dp = sqrt(var(veloc)), vari = var(veloc), 
               cv = 100*((sqrt(var(veloc))/mean(veloc))), 
               ca = skewness(veloc), minimo = min(veloc), maximo = max(veloc))
xtable(cturb)


# boxplot
boxplot(veloc~grupo_fac, names = c("1","2","3","4","5"), xlab = "Tipo de turbina", 
        ylab = "Tempo de vida (em milhões de ciclos)", cex = 1.2, cex.axis = 1.2, 
        cex.lab = 1.2)



#####################
#####################
# Análise Inferencial
#####################
#####################


################
################
# Modelo Normal

fit.model <- lm(veloc~grupo_fac)
summary(fit.model)
anova(fit.model)
AIC(fit.model)
BIC(fit.model)
X11()
diagnorm(fit.model)
envelnorm(fit.model)


########################
########################
# Modelo Gama (link log)

fit.model <- glm(veloc~grupo_fac, family = Gamma(link = "log"))
summary(fit.model)


# Usando a estimativa de MV de phi:
summary(fit.model, dispersion = gamma.dispersion(fit.model))
xtable(summary(fit.model, dispersion = gamma.dispersion(fit.model)))


# Estimativa do parâmetro phi (é o inverso da dispersão):

## Default R (por padrão, a estimativa de phi é dada pelo Método dos Momentos)
ephi <- 1/summary(fit.model)$dispersion # valor estimado pelo default do glm

## Método dos momentos
ro <- resid(fit.model, type = "response") # resíduos do modelo
ephiM <-(45)/sum((ro/(fitted(fit.model)))^2) # = ephi
                                             # 45 = n - p, onde n = 50 e p (numero de paramtros) = 5 (alpha, beta_2, beta_3, beta_4, beta_5)
                                             # sum((ro/(fitted(fit.model)))^2 = soma dos resíduos ao quadrado dividido por Mi, que seria o Vi, pra gente fazer a soma do desvio ao quadrado

## MV
resultphiMV <- gamma.shape(fit.model) # resultado da estimativa via MV
ephiMV <- resultphiMV$alpha # valor da estimativa de phi por MV 
sephiMV <- resultphiMV$SE # valor do desvio padrão por MV

ez <- qnorm(0.975)
ICphi <- cbind(ephiMV-ez*sephiMV,ephiMV+ez*sephiMV) # intervalo de confiança assintótico da estimativa de phi

## Outro exemplo de estimador de phi (Não consistente): (n-p)/função desvio do modelo
ephiNC <- 45/deviance(fit.model)


# AIC e BIC (outra forma de extrair os valores de AIC e BIC do modelo, com base no phi de MV)
AICBICgama(veloc, fitted(fit.model), ephiMV, 50, 6)


# Deviance do modelo
desv <- deviance(fit.model)*ephiMV # função desvio escalonada, pois esta multiplicando por phi de MV
pvalordesv <- 1-pchisq(desv, df = 45) # aqui estamos utilizando a dist. qui-quadrado como referência ; lembrando que isso vale pra quando
                                      # o valor de phi é suficientemente grande, para aproximar para uma dist. qui-quadrado
                                      # p-valor = 0.233, indicando que o modelo se ajustou de modo satisfatório aos dados


# Gráficos de diagnóstico
diaggama(fit.model, 1) # em comparação com o modelo Normal, o ponto discrepante está mais proximo, ele parece ter se adequado mais
                       # em relação ao histograma, ele ta mais simétrico em torno do zero
                       # outro gráfico interessante de avaliar no MLG, é o grafico do "predito linear vs variavel z" ; esse gráfico é pra identificar se a função de ligação
                       # está adequada ao modelo ; o ideal é que esse grafico fique proximo de um comportamento/estrutura linear
envelgama(fit.model, "log", 1) # aquele ponto discrepante que ficava fora das bandas, agr ja esta dentro ; apresenta um comportamento melhor que o modelo anterior


# Medias preditas pelo modelo
resmupred <- predict(fit.model, type = c("response"), se.fit = TRUE) # valores preditos da media e do desvio padrao para cada turbina
mupred1 <- unique(resmupred$fit) # valores preditos da media de cada turbina
epmupred1 <- unique(resmupred$se.fit) # valores preditos do desvio padrao de cada turbina



###############################
###############################
# Modelo Gama (link identidade)

fit.model <- glm(veloc~grupo_fac, family = Gamma(link = "identity"))
summary(fit.model)


# Usando a estimativa de MV de phi
summary(fit.model, dispersion = gamma.dispersion(fit.model))
xtable(summary(fit.model, dispersion = gamma.dispersion(fit.model)))


# Estimativa do parãmetro phi

## MV
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 
sephiMV <- resultphiMV$SE

ez <- qnorm(0.975)
ICphi <- cbind(ephiMV-ez*sephiMV,ephiMV+ez*sephiMV)

# Estimativas dos parâmetros
resultmg2 <- coef(summary(fit.model, dispersion = gamma.dispersion(fit.model)))

ez <- qnorm(0.975)
ICbeta <-cbind(resultmg2[,1]-ez*resultmg2[,2], resultmg2[,1]+ez*resultmg2[,2])

resultmg2 <- cbind(resultmg2[,1], resultmg2[,2], ICbeta,resultmg2[,3], resultmg2[,4])
resultmg2 <- rbind(resultmg2, cbind(c(ephiMV), c(sephiMV), ICphi, 0, 0))
xtable(resultmg2)


# AIC e BIC
AICBICgama(veloc, fitted(fit.model), ephiMV, 50, 6)


# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=45)


# Gráficos de diagnóstico
X11()
diaggama(fit.model, 1) # graficos bem parecidos com o q a gente viu na funçao logaritmica
envelgama(fit.model, "identity", 1) # aparenta estar melhor, pois agora todos os pontos estao dentro das bandas
                                    # a função identidade aparentou se adequar melhor aos dados

# Medias preditas pelo modelo
resmupred <- predict(fit.model,type = c("response"),se.fit=TRUE)
mupred2 <- unique(resmupred$fit)
epmupred2 <- unique(resmupred$se.fit)

ez <-qnorm(0.975)
plotCI(c(1,2,3,4,5), mupred1, uiw = epmupred1*ez, liw = epmupred1*ez, pch = 19,lwd = 2,
       col = 1, xlab = "tipo de turbina", ylab = "tempo de vida médio", cex = 1.2,
       cex.lab = 1.2, cex.axis = 1.2, xlim = c(1,5.5))
plotCI(c(1.2,2.2,3.2,4.2,5.2), mupred2, uiw = epmupred2*ez, liw = epmupred2*ez, pch = 17,
       cex = 1.2, lwd = 2, col = 2, add = TRUE)
lines(c(0.9,1.9,2.9,3.9,4.9), cturb$media, pch = 15, lwd = 2, col = 3, type = "p", cex = 1.2)
legend(1.5, 18, c("observada","modelo gama - ligação log","modelo gama - ligação identidade"), 
       col = c(3,1,2), pch = c(15,19,17), bty = "n", cex = 1.2)



###################################
###################################
# Modelo reduzido (link identidade)


#grupo_facr <- revalue(grupo_fac,c("1"="1","2"="2","3"="1","4"="1","5"="1"))
grupo_facr <- revalue(grupofacr, c("1"="1","2"="2","3"="1","4"="1","5"="5"))
fit.model <- glm(veloc~grupo_facr, family = Gamma(link = "identity"))
summary(fit.model)


# Usando a estimativa de MV de phi
summary(fit.model, dispersion = gamma.dispersion(fit.model))
xtable(summary(fit.model, dispersion = gamma.dispersion(fit.model)))


# Estimativa do parâmetro phi

## MV
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 
sephiMV <- resultphiMV$SE

ez <- qnorm(0.975)
ICphi <- cbind(ephiMV-ez*sephiMV, ephiMV+ez*sephiMV)

## Estimativas dos parâmetros ( os betas)
resultmg3 <- coef(summary(fit.model, dispersion = gamma.dispersion(fit.model)))

ez <- qnorm(0.975)
ICbeta <-cbind(resultmg3[,1]-ez*resultmg3[,2], resultmg3[,1]+ez*resultmg3[,2])

resultmg3 <- cbind(resultmg3[,1], resultmg3[,2], ICbeta, resultmg3[,3], resultmg3[,4])
resultmg3 <- rbind(resultmg3, cbind(c(ephiMV), c(sephiMV), ICphi, 0, 0)) # a ultima linha se refere ao parâmetro de precisão phi
xtable(resultmg3)


# AIC e BIC
AICBICgama(veloc, fitted(fit.model), ephiMV, 50, 4) # possui valores menores, ou seja, apresentou melhora


# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=50-3)


# Gráficos de diagnóstico
X11()
diaggama(fit.model, 1) # no histograma ja é possivel notar um comportamento mais simetrico
envelgama(fit.model, "identity", 1)


# Médias preditas pelo modelo
resmupred <- predict(fit.model, type = c("response"), se.fit = TRUE)
mupred3 <- unique(resmupred$fit)
mupred3 <- c(mupred3[1], mupred3[2], mupred3[1], mupred3[1], mupred3[3])
epmupred3 <- unique(resmupred$se.fit)
epmupred3 <- c(epmupred3[1], epmupred3[2], epmupred3[1], epmupred3[1], epmupred3[3])

ez <-qnorm(0.975)
plotCI(c(1,2,3,4,5), mupred3, uiw = epmupred3*ez, liw = epmupred3*ez, pch = 19, lwd = 2,
       col = 1, xlab = "tipo de turbina", ylab = "tempo de vida médio", cex = 1.2,
       cex.lab = 1.2, cex.axis = 1.2, xlim = c(1,5.5))
lines(c(0.9,1.9,2.9,3.9,4.9), cturb$media, pch = 15, lwd = 2, col = 3, type = "p", cex = 1.2)
legend(1.5, 18, c("observada","modelo gama 3"), col = c(3,1), pch = c(15,19), bty = "n", cex = 1.2)
xtable(cbind(mupred3,epmupred3,mupred3-ez*epmupred3,mupred3+ez*epmupred3))


# Análise preditiva
par(mfrow=c(1,1))
boxplot(veloc~grupo_fac, names = c("1","2","3","4","5"), xlab = "tipo de turbina", 
        ylab = "tempo de vida (em milhões de ciclos)", cex = 1.2, cex.axis = 1.2, cex.lab = 1.2,
        ylim = c(0,40), outwex = 2)

npred <- 1000
velocpred <- c((mupred3[1]/ephiMV)*rgamma(npred,ephiMV))

for (j in 2:5)
{
  velocpred <- c(velocpred, c((mupred3[j]/ephiMV)*rgamma(npred,ephiMV)))
}

auxgrupo <- c(rep(1,npred), rep(2,npred), rep(3,npred), rep(4,npred), rep(5,npred))
auxgrupofac <- factor(auxgrupo, c("1","2","3","4","5"))
boxplot(velocpred~auxgrupofac, names = c("1","2","3","4","5"), xlab = "tipo de turbina",
        ylab = "tempo de vida (em milhões de ciclos)", cex = 1.2, cex.axis = 1.2, cex.lab = 1.2,
        add = TRUE, border = 2)
legend(2, 37, c("observada","predita"), col = c(1,2), lty = c(1,1), lwd = c(2,2), bty = "n", cex = 1.2)