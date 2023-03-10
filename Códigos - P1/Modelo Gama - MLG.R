####################### MLG - MODELO GAMA ###########################

# Definindo variaveis:

n <- 100 # tamanho da amostra
x <- runif(n,1,10) # gerando valores uniformes para a variavel explicativa
eta <- 2 - 0.6*x # preditor linear com beta_0 = 2 e beta_1 = -0.6

mu <- exp(eta) # media, ja que estamos usando a fun??o log, log(mi) = eta, ent?o, mi = exp(eta)
phi <- 2.5 # parametro verdadeiro
k <- mu/phi
r <- phi

y <- NULL
for(i in 1:n){ # gerando dados da fun??o gama
  y[i]=rgamma(1,shape = r,scale = k[i])}

X <- cbind(rep(1,n),x) # matriz de delineamento com intercepto (beta_0), e a propria variavel explicativa x


# Agora, vamos aplicar o Algoritmo do Escore de Fisher (AEF):

N <- 10 # quantidade de itera??es
fi <- numeric() # estimativa de phi
beta <- matrix(0, nrow=N, ncol=2) # estimativa dos coeficientes do modelo linear (beta_0 e beta_1)
                                  # nrow = 10, pois beta vai guardar as estimativas de b0 e b1 a cada itera??o
epsilon <- numeric() # epsilon vai guardar a medida de diferen?a quadr?tica entre estimativas obtidas para 
                     # duas itera??es sucessivas, ou seja, vai calcular a distancia entre o beta_i e o beta_i-1

## Vamos realizar 10 itera??es:

for(i in 1:N){
  if(i==1) { # Se for a primeira itera??o, ent?o usamos nossos chutes iniciais
    mu <- rep(mean(y),n) 
    eta <- log(mu)
    fi <- mean(y)^2/var(y)}
  
  if (i!=1) {
    eta <- X %*% beta[i-1,] # Preditor linear calculado com as estimativas do passo anterior
    mu <- exp(eta) # a partir da segunda itera??o usamos os resultados do passo anterior para obter estimativas no passo atual
  }
  
  # calculo dos betas estimados por MV ; atualia??o dos betas:
  vmu <- mu^2 # fun??o de vari?ncia (V_i) avaliada em mu
  glinhamu <- 1/(mu) # derivada da fun??o de liga??o avaliadas em mu
  
  z <- eta+(y-mu)*glinhamu # Vetor z avaliada em mu
  W <- diag(1,n) # matriz identidade ; matriz diagonal avaliada em mu
  beta[i,] <- solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%z) # Solu??o de minimos quadrados ponderados para o vetor beta no passo i.
  
  # calculo do fi estimado por MV ; atualiza??o do fi:
  theta <- -1/mu
  b <- -log(-theta)
  clinha <- sum(log(y))+n*(log(fi)+1-digamma(fi))
  c2linha <- n*(1/fi-trigamma(fi))
  fi <- fi-(sum(y*theta-b)+clinha)/(c2linha)
  
  if(i>1) epsilon[i-1]=sum(((beta[i,]-beta[i-1,])/beta[i-1,])**2) # somatorio da distancia relativa ao quadrado
  
}

## OBS: normalmente, colocamos um criterio de parada no if feito acima ; e esse criterio de parada ? o epsilon
#       se ele for muito pequeno ( < 0.0001, por exemplo), significa que a convergencia foi atingida, e n?o 
#       ? necessario realizar mais itera??es


# Plotando o ajuste:

plot(x, y, pch = 20, cex = 1.5, las = 1)

ajuste = function(x) exp(beta[10,1]+beta[10,2]*x) # exp(eta)
curve(ajuste, 0, 40, add = T, col =2 , lwd = 2)


# Vamos tentar o ajuste de uma regress?o linear simples:

ajuste2 <- lm(y ~ x)
abline(coefficients(ajuste2), col = 'blue', lwd = 2) # nota-se que n?o consegue achar essa curva exponencial
                                                     # em algumas partes os valores est?o subestimados, em outras os valores est?o superestimados
                                                     # n?o se ajusta adequadamente


# Agora, vamos ajustar o modelo declarando a log-verossimilhan?a a um otimizador do R:

require(bbmle)

logvero = function(b0, b1, phi)
  -sum(dgamma(y, shape = phi, scale = exp(b0+b1*x)/phi, log = T)) # logvero armazena a fun??o de log verossimilhan?a (-)

ajuste <- lm(log(y)~x) # Regressando g(y/n) em fun??o de x ; ajuste do modelo de regressao para obter valores iniciais para b0 e b1

est2 <- mle2(logvero, start = list(b0 = ajuste$coefficients[1],
                                   b1 = ajuste$coefficients[2], 
                                   phi = 1/summary(ajuste)$sigma))
est2


# Finalmente, usemos a fun??o glm (fun??o do R para MLG):

ajuste <- glm(y~x, family = Gamma(link = "log"))
names(ajuste)
coef(ajuste) # Estimativas dos coeficientes do modelo
fitted(ajuste) # M?dias ajustadas pelo modelo para cada valor de x na amostra
predict(ajuste, newdata = data.frame(x = c(3.5,5.5,7.5))) # Estimativas (na escala do preditor) para x=3.5; x=5.5 e x=7.5
predict(ajuste, newdata = data.frame(x = c(3.5,5.5,7.5)), type = 'response') # Estimativas (na escala da resposta - probabilidades estimadas) para x=3.5; x=5.5 e x=7.5
summary(ajuste) # Resumo do modelo ajustado contendo, dentre outras coisas, as estimativas dos betas e os correspondnetes erros padr?es
vcov(ajuste) # Matriz de vari?ncias e covari?ncias estimada

## OBS: o glm fornece o estimador de momentos
#       para obter o estimador de MV, basta rodar o codigo, para esse caso, gamma.shape(ajuste), do pacote MASS