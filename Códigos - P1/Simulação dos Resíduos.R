# Dados heteroced?sticos
source("https://www.ime.unicamp.br/~cnaber/diag_gama.r")
source("https://www.ime.unicamp.br/~cnaber/envel_gama.r")

# Fun??es auxiliares para o modelo gama
source("https://www.ime.unicamp.br/~cnaber/AuxModGama.r")

# pacotes
library(MASS)
library(sn)



n <- 200
x <- runif(n,0,1)
#aux <- c(x[1:50],1/x[51:100])
beta0 <- 1
beta1 <- 1.2
pred <- beta0+beta1*x # preditor linear
phi <- 5



# Modelo 1
# Igual: modelo simulado-modelo ajustado
y <- rgamma(n, phi)
mu <- exp(pred)
y <- (mu/phi)*y

fit.model <- glm(y~x, family = Gamma(link = "log"))
diaggama(fit.model, 1) # valores apresenta comportamento assim?trico em torno do zero com varia??o constantes
                       # com alguns pontos fora da banda, mas sem ultrapassar o limite de 5%
                       # grafico inferior direito apresenta comportamento linear, apresentar uma estrutura linear
envelgama(fit.model, "log", 1) # comportamento adequado com todos os pontos dentro da banda

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 2
# Diferentes phi's
vphi <- exp(3*x)
y <- rgamma(n,vphi)
mu <- exp(pred)
y <- (mu/vphi)*y

fit.model <- glm(y~x, family = Gamma(link = "log"))
diaggama(fit.model,1) # no grafico superior direito nota-se uma heterocedasticidade que o residuo n?o captou ; mo come?o tem-se uma dispers?o alta q vai diminuindo
                      # a msm analisa feita acima tambem vale para o grafico inferior direito
                      # em rela??o ao grafico inferior esquerdo, nota-se que ele est? mais concentrado, com uma menor dispers?o 
envelgama(fit.model,"log",1) # ja nota-se pontos saindo da banda na parte sperior e inferior

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV < -gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 3
# Fun??o de liga??o log e ajuste liga??o identidade
y <- rgamma(n, phi)
mu <- exp(pred)
y <- (mu/phi)*y

fit.model <- glm(y~x, family = Gamma(link = "identity"))
diaggama(fit.model, 1) # os graficos de cima similar ao modelo 1, comportamento aleatorio em torno do zero com variancia constante, porem com mais pontos discrepantes
                       # no grafico do preditor linear a gente tamb?m observa alguns pontos mais distantes da reta
envelgama(fit.model, "identity", 1) # nota-se uma leve sa?da das bandas

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 4
# Outra distribui??o (Normal assim?trica)
beta0 <- 3
beta1 <- 1.1
pred <- beta0+beta1*x # preditor linear
phi <- 5
mu <- exp(pred) # segue a mesma estrutura da fun??o de liga??o, da estrutura do modelo
y <- rsn(n,mu,phi,-20)
hist(y)

fit.model <- glm(y~x, family = Gamma(link = "log"))
diaggama(fit.model, 1) # em rela??o ao gr?fico do preditor linear, aparenta ter comportamento linear
                       # em rela??o ao grafico superior esquerdo, nota-se um comportamento aleatorio em torno do zero com variancia constante, porem com alguns pontos discrepantes
                       # em rela??o ao grafico superior direito, nota-se uma certa heterocedasticidade q n foi captada pelo modelo, come?a com uma variancia alta e esta vai diminuindo
                       # em rela??o ao gr?fico inferior esquerdo, nota-se um comportamento bastante assimetrico
envelgama(fit.model, "log", 1) # nota-se q muitos pontos estao fora da banda

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV<-gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 5
# Outra distribui??o (t de student assim?trica)
beta0 <- 3
beta1 <- 1.1
pred <- beta0+beta1*x
phi <- 5
mu <- exp(pred)
y <- rsn(n,mu,phi,-20,3)
hist(y)

fit.model <- glm(y~x, family = Gamma(link = "log"))
diaggama(fit.model, 1) # em rela??o ao gr?fico do preditor linear, aparenta ter comportamento linear, porem tamb?m e percept?vel uma heterocedasticidade
                       # em rela??o ao grafico superior esquerdo, nota-se um comportamento aleatorio em torno do zero com variancia constante, porem com alguns pontos discrepantes
                       # em rela??o ao grafico superior direito, nota-se uma certa heterocedasticidade q n foi captada pelo modelo, come?a com uma variancia alta e esta vai diminuindo
                       # em rela??o ao gr?fico inferior esquerdo, nota-se um comportamento assimetrico
envelgama(fit.model, "log", 1) # tbm apresenta pontos fora da banda, porem menos que o modelo 4

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)