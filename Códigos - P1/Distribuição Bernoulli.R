### Exemplo - Distribui??o de Bernoulli ###


y <- c(0,0,0,0,1,0,1,1,1,1) # onde 0 = fracasso e 1 = sucesso
x <- 1:10
plot(x,y,pch=20) # ou seja, a medida de que x aumenta, maior a probabilidade de sucesso

X <- matrix(c(rep(1,10),x), nrow = 10) # matriz do modelo


# Implementando o Algor?tmo Escore de Fisher (AEF): oi

beta <- matrix(0, nrow = 10, ncol = 2) # beta vai guardar as estimativas de b0 e b1 a cada itera??o
epsilon <- numeric() # epsilon vai guardar a medida de diferen?a quadr?tica entre as estimativas obtidas para duas itera??es sucessivas

## Vamos realizar 10 itera??es

for(i in 1:10){
  if(i==1) {mu <- c(0.1,0.1,0.1,0.1,0.9,0.1,0.9,0.9,0.9,0.9)
  ### Se for a primeira itera??o, ent?o usamos nossos chutes iniciais
  eta <- log(mu/(1-mu))} ### preditor linear
  
  if (i!=1) { ### atualiza??o do eta
    eta <- X %*% beta[i-1,] ### Preditor linear calculado com as estimativas do passo anterior
    mu <- exp(eta)/(exp(eta)+1) ### a partir da segunda itera??o, usamos os resultados do passo anterior para obter estimativas no passo atual
  }
  
  ### Fun??o de variancia e derivada da  fun??o de liga??o avaliadas em mu
  vmu <- mu*(1-mu)
  glinhamu <- 1/(mu*(1-mu))
  
  ### Vetor z e matriz diagonal avaliadas em mu
  z <- eta+(y-mu)*glinhamu
  W <- diag(as.numeric((vmu*(glinhamu**2))**(-1)))
  
  beta[i,] <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% z) ### Solu??o de minimos quadrados ponderados para o vetor beta no passo i
  if(i>1) epsilon[i-1]=sum(((beta[i,]-beta[i-1,])/beta[i-1,])**2)
  
}


# Plotando o ajuste.

plot(x,y,pch=20,cex=1.5,ylim=c(-0.2,1.2),las=1)
ajuste=function(x) exp(beta[10,1]+beta[10,2]*x)/(1+exp(beta[10,1]+beta[10,2]*x))
curve(ajuste,0,40,add=T)


# Vamos tentar o ajuste de uma regress?o linear simples:
ajuste2 <- lm(y ~ x)
abline(coefficients(ajuste2), col='red')
## OBS: FICA RUIM, COME?A A TER VALORES PREDITOS FORA DO INTRVALO (0,1)


# Agora, vamos ajustar o modelo declarando a log-verossimilhan?a a um otimizador do R:

require(bbmle)
logvero=function(b0,b1)
  -sum(dbinom(y,1,exp(b0+b1*x)/(1+exp(b0+b1*x)), log=T))
## OBS: logvero armazena a fun??o de log verossimilhan?a (-)

ylinha <- c(0.1,0.1,0.1,0.1,0.9,0.1,0.9,0.9,0.9,0.9)
ajuste <- lm(log((ylinha)/(1-ylinha))~x) ### Regressando g(y/n) em fun??o de x para obter valores iniciais para b0 e b1
est2 <- mle2(logvero,start=list(b0=ajuste$coefficients[1],b1=ajuste$coefficients[2]))
est2 # os valores est?o muito proximos da estima??o do AEF


# Finalmente, usemos a fun??o glm.

ajuste <- glm(y~x,family=binomial(link = "logit"))
names(ajuste)
coef(ajuste) ### Estimativas dos coeficientes do modelo ; ? o mesmo resultado do AEF
fitted(ajuste) ### Probabilidades ajustadas pelo modelo para cada valor de x na amostra
predict(ajuste,newdata=data.frame(x=c(3.5,5.5,7.5))) ### Estimativas (na escala do preditor) para x=3.5; x=5.5 e x=7.5
predict(ajuste,newdata=data.frame(x=c(3.5,5.5,7.5)),type='response') ### Estimativas (na escala da resposta - probabilidades estimadas) para x=3.5; x=5.5 e x=7.5
summary(ajuste) ### Resumo do modelo ajustado contendo, dentre outras coisas, as estimativas dos betas e os correspondnetes erros padr?es
vcov(ajuste) ### Matriz de vari?ncias e covari?ncias estimada

sqrt(vcov(ajuste)[1,1]) # erro padr?o
sqrt(vcov(ajuste)[2,2]) # erro padr?o
ajuste$coefficients[1]/sqrt(vcov(ajuste)[1,1]) # estat?stica z
2*(1 - pnorm(ajuste$coefficients[1]/sqrt(vcov(ajuste)[1,1]))) # p-valor


######################################################################################


### Outro Exemplo - Dist. Bernoulli ###

n <- 100
x <- runif(n,1,10)
eta <- 2-.6*x # preditor linear, onde intercepto = 2 e coeficiente angular = -0.6
mu <- exp(eta)/(1+exp(eta)) # fun??o de liga??o = inversa da fun??o logito ; probabilidade de sucesso
y <- rbinom(n,1,mu)
plot(x, y, pch=20)

X <- matrix(c(rep(1,n),x),nrow=n) ### Matriz do modelo


# Implementando o Algor?tmo Escore de Fisher (AEF):

N <- 10

beta <- matrix(0, nrow=N, ncol=2)
epsilon <- numeric()
### beta vai guardar as estimativas de b0 e b1 a cada itera??o
### Epsilon vai guardar a medida de diferen?a quadr?tica entre estimativas
### obtidas para duas itera??es sucessivas


# Vamos realizar 10 itera??es

for(i in 1:N){
  if(i==1) {mu<-rep(.1,n);mu[y==1]<-.9
  ### Se for a primeira itera??o ent?o usamos nossos chutes iniciais
  eta <- log(mu/(1-mu))}
  
  if (i!=1) {
    eta <- X %*% beta[i-1,] ### Preditor linear calculado com as estimativas do passo anterior
    mu <- exp(eta)/(exp(eta)+1) ### a partir da segunda itera??o, usamos os resultados do passo anterior para obter estimativas no passo atual
  }
  
  vmu <- mu*(1-mu)
  glinhamu <- 1/(mu*(1-mu))
  ### Fun??o de variancia e derivada da fun??o de liga??o avaliadas em mu
  
  z <- eta+(y-mu)*glinhamu
  W <- diag(as.numeric((vmu*(glinhamu**2))**(-1)))
  ### Vetor z e matriz diagonal avaliadas em mu
  
  beta[i,] <- solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%z) ### Solu??o de minimos quadrados ponderados para o vetor beta no passo i
  if(i>1) epsilon[i-1]=sum(((beta[i,]-beta[i-1,])/beta[i-1,])**2)
  
}


# Plotando o ajuste.

plot(x,y,pch=20,cex=1.5,ylim=c(-0.2,1.2),las=1)
ajuste=function(x) exp(beta[N,1]+beta[N,2]*x)/(1+exp(beta[N,1]+beta[N,2]*x))
curve(ajuste,0,40,add=T)


# Vamos tentar o ajuste de uma regress?o linear simples:
ajuste2 <- lm(y ~ x)
abline(coefficients(ajuste2), col='red')
## OBS: FICOU HORRIVEL

# Agora, vamos ajustar o modelo declarando a log-verossimilhan?a a um otimizador do R:

require(bbmle)
logvero=function(b0,b1)
  -sum(dbinom(y,1,exp(b0+b1*x)/(1+exp(b0+b1*x)), log=T))
## logvero armazena a fun??o de log verossimilhan?a (-)

ylinha <- rep(.1,n); ylinha[y==1]<-.9
ajuste <- lm(log((ylinha)/(1-ylinha))~x) ### Regressando g(y/n) em fun??o de x para obter valores iniciais para b0 e b1
est2 <- mle2(logvero,start=list(b0=ajuste$coefficients[1],b1=ajuste$coefficients[2]))
est2 # estimativas muito proximas do AEF


# Finalmente, usemos a fun??o glm.

ajuste <- glm(y~x,family=binomial(link = "logit"))
names(ajuste)
coef(ajuste) ### Estimativas dos coeficientes do modelo ; os mesmos valores que o do AEF
fitted(ajuste) ### Probabilidades ajustadas pelo modelo para cada valor de x na amostra
predict(ajuste,newdata=data.frame(x=c(3.5,5.5,7.5))) ### Estimativas (na escala do preditor) para x=3.5; x=5.5 e x=7.5
predict(ajuste,newdata=data.frame(x=c(3.5,5.5,7.5)),type='response') ### Estimativas (na escala da resposta - probabilidades estimadas) para x=3.5; x=5.5 e x=7.5
summary(ajuste) ### Resumo do modelo ajustado contendo, dentre outras coisas, as estimativas dos betas e os correspondnetes erros padr?es
vcov(ajuste) ### Matriz de vari?ncias e covari?ncias estimada

sqrt(vcov(ajuste)[1,1]) # erro padr?o
sqrt(vcov(ajuste)[2,2]) # erro padr?o
ajuste$coefficients[1]/sqrt(vcov(ajuste)[1,1]) # estat?stica z
2*(1 - pnorm(ajuste$coefficients[1]/sqrt(vcov(ajuste)[1,1]))) # p-valor