#------------------------------------------------------------#
# Para  rodar este programa deixe no objeto fit.model a sa?da 
# do ajuste da regress?o do modelo normal linear. Deixe tamb?m
# os dados dispon?veis atrav?s do comando attach(...). Depois
# use o comando source(...) no R ou S-plus para executar o 
# programa. A sequ?ncia de comandos ? a seguinte:
#
#       fit.model <- ajuste
#       attach(dados)
#       source("envel_norm")
#
# A sa?da ser? o gr?fico de envelope para o res?duo
# padronizado. Para colocar  um  t?tulo no gr?fico ap?s a
# sa?da use o comando title("...").
#------------------------------------------------------------#
par(mfrow=c(1,1))

X <- model.matrix(fit.model) # fit.model = modelo ajustado com a fun??o lm ; x = matriz de planejamento, ou seja, das vari?veis preditoras
n <- nrow(X) # tamanho da amostra ; numero de linhas da matriz
p <- ncol(X) # quantidade de parametros ; numero de colunas da matriz
H <- X%*%solve(t(X)%*%X)%*%t(X) # matriz H
h <- diag(H) # h = elementos da diagonal da matriz H
si <- lm.influence(fit.model)$sigma # si = sigma ao quadrado chapeu com a informa??o i deletada
r <- resid(fit.model) # residuo padronizado
tsi <- r/(si*sqrt(1-h)) # residuo studentizado

ident <- diag(n) # matriz identidade

# AQUI COME?A O PROCESSO DE SIMULA??O
epsilon <- matrix(0,n,100)
e <- matrix(0,n,100) # constroi a matriz com m = 100, ou seja, 100 repeti??es, e n ? o tamanho da amosra ; entao temos uma matriz com n linhas e 100 colunas
e1 <- numeric(n) # vetor de limite de tamanho n
e2 <- numeric(n) # vetor de limite de tamanho n

##
for(i in 1:100){ # aqui estamos varrendo as colunas
  epsilon[,i] <- rnorm(n,0,1) # aqui geramos a n observa??es normal padrao
  e[,i] <- (ident - H)%*%epsilon[,i] # aqui geramos o r* (segundo passo da cria??o do grafico envelope)
  u <- diag(ident - H)
  e[,i] <- e[,i]/sqrt(u) # aqui geramos o t*_i, que o res?duo padronizado
  e[,i] <- sort(e[,i]) } # aqui ordemos os t*_i's

##
for(i in 1:n){ # aqui estamos varrendo as linhas, ou seja, pra cada valor na amostra
  eo <- sort(e[i,]) # aqui se trata do residuo ordenado
  e1[i] <- (eo[2]+eo[3])/2 # aqui ? o calculo do limite inferior ; pega-se as posi??es 2 e 3, para se obter o percentil de 2.5
  e2[i] <- (eo[97]+eo[98])/2 } # aqui ? o calculo do limite superior ; pega-se as posi??es 97 e 98, para se obter o percentil de 97.5

##
med <- apply(e,1,mean) # linha de refer?ncia (setimo passo da cria??o do grafico envelope)
faixa <- range(tsi,e1,e2)

##
par(pty="s")
qqnorm(tsi,xlab="Percentil da N(0,1)",
       ylab="Residuo Studentizado", ylim=faixa, pch=16, main="")
par(new=TRUE)
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1, main="")
par(new=TRUE)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1, main="")
par(new=TRUE)
qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2, main="")
#------------------------------------------------------------#