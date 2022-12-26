# Dados heterocedásticos
source("https://www.ime.unicamp.br/~cnaber/diag_gama.r")
source("https://www.ime.unicamp.br/~cnaber/envel_gama.r")

# Funções auxiliares para o modelo gama
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
diaggama(fit.model, 1) # valores apresenta comportamento assimétrico em torno do zero com variação constantes
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
diaggama(fit.model,1) # no grafico superior direito nota-se uma heterocedasticidade que o residuo não captou ; mo começo tem-se uma dispersão alta q vai diminuindo
                      # a msm analisa feita acima tambem vale para o grafico inferior direito
                      # em relação ao grafico inferior esquerdo, nota-se que ele está mais concentrado, com uma menor dispersão 
envelgama(fit.model,"log",1) # ja nota-se pontos saindo da banda na parte sperior e inferior

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV < -gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 3
# Função de ligação log e ajuste ligação identidade
y <- rgamma(n, phi)
mu <- exp(pred)
y <- (mu/phi)*y

fit.model <- glm(y~x, family = Gamma(link = "identity"))
diaggama(fit.model, 1) # os graficos de cima similar ao modelo 1, comportamento aleatorio em torno do zero com variancia constante, porem com mais pontos discrepantes
                       # no grafico do preditor linear a gente também observa alguns pontos mais distantes da reta
envelgama(fit.model, "identity", 1) # nota-se uma leve saída das bandas

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 4
# Outra distribuição (Normal assimétrica)
beta0 <- 3
beta1 <- 1.1
pred <- beta0+beta1*x # preditor linear
phi <- 5
mu <- exp(pred) # segue a mesma estrutura da função de ligação, da estrutura do modelo
y <- rsn(n,mu,phi,-20)
hist(y)

fit.model <- glm(y~x, family = Gamma(link = "log"))
diaggama(fit.model, 1) # em relação ao gráfico do preditor linear, aparenta ter comportamento linear
                       # em relação ao grafico superior esquerdo, nota-se um comportamento aleatorio em torno do zero com variancia constante, porem com alguns pontos discrepantes
                       # em relação ao grafico superior direito, nota-se uma certa heterocedasticidade q n foi captada pelo modelo, começa com uma variancia alta e esta vai diminuindo
                       # em relação ao gráfico inferior esquerdo, nota-se um comportamento bastante assimetrico
envelgama(fit.model, "log", 1) # nota-se q muitos pontos estao fora da banda

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV<-gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)



# Modelo 5
# Outra distribuição (t de student assimétrica)
beta0 <- 3
beta1 <- 1.1
pred <- beta0+beta1*x
phi <- 5
mu <- exp(pred)
y <- rsn(n,mu,phi,-20,3)
hist(y)

fit.model <- glm(y~x, family = Gamma(link = "log"))
diaggama(fit.model, 1) # em relação ao gráfico do preditor linear, aparenta ter comportamento linear, porem também e perceptível uma heterocedasticidade
                       # em relação ao grafico superior esquerdo, nota-se um comportamento aleatorio em torno do zero com variancia constante, porem com alguns pontos discrepantes
                       # em relação ao grafico superior direito, nota-se uma certa heterocedasticidade q n foi captada pelo modelo, começa com uma variancia alta e esta vai diminuindo
                       # em relação ao gráfico inferior esquerdo, nota-se um comportamento assimetrico
envelgama(fit.model, "log", 1) # tbm apresenta pontos fora da banda, porem menos que o modelo 4

# Estimativa de phi po MV e seu respectivo erro padrao
resultphiMV <- gamma.shape(fit.model)
ephiMV <- resultphiMV$alpha 

# Deviance do modelo
desv <- deviance(fit.model)*ephiMV
pvalordesv <- 1-pchisq(desv,df=198)