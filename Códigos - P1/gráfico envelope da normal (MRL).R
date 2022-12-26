#------------------------------------------------------------#
# Para  rodar este programa deixe no objeto fit.model a saída 
# do ajuste da regressão do modelo normal linear. Deixe também
# os dados disponíveis através do comando attach(...). Depois
# use o comando source(...) no R ou S-plus para executar o 
# programa. A sequência de comandos é a seguinte:
#
#       fit.model <- ajuste
#       attach(dados)
#       source("envel_norm")
#
# A saída será o gráfico de envelope para o resíduo
# padronizado. Para colocar  um  título no gráfico após a
# saída use o comando title("...").
#------------------------------------------------------------#
par(mfrow=c(1,1))

X <- model.matrix(fit.model) # fit.model = modelo ajustado com a função lm ; x = matriz de planejamento, ou seja, das variáveis preditoras
n <- nrow(X) # tamanho da amostra ; numero de linhas da matriz
p <- ncol(X) # quantidade de parametros ; numero de colunas da matriz
H <- X%*%solve(t(X)%*%X)%*%t(X) # matriz H
h <- diag(H) # h = elementos da diagonal da matriz H
si <- lm.influence(fit.model)$sigma # si = sigma ao quadrado chapeu com a informação i deletada
r <- resid(fit.model) # residuo padronizado
tsi <- r/(si*sqrt(1-h)) # residuo studentizado

ident <- diag(n) # matriz identidade

# AQUI COMEÇA O PROCESSO DE SIMULAÇÃO
epsilon <- matrix(0,n,100)
e <- matrix(0,n,100) # constroi a matriz com m = 100, ou seja, 100 repetições, e n é o tamanho da amosra ; entao temos uma matriz com n linhas e 100 colunas
e1 <- numeric(n) # vetor de limite de tamanho n
e2 <- numeric(n) # vetor de limite de tamanho n

##
for(i in 1:100){ # aqui estamos varrendo as colunas
  epsilon[,i] <- rnorm(n,0,1) # aqui geramos a n observações normal padrao
  e[,i] <- (ident - H)%*%epsilon[,i] # aqui geramos o r* (segundo passo da criação do grafico envelope)
  u <- diag(ident - H)
  e[,i] <- e[,i]/sqrt(u) # aqui geramos o t*_i, que o resíduo padronizado
  e[,i] <- sort(e[,i]) } # aqui ordemos os t*_i's

##
for(i in 1:n){ # aqui estamos varrendo as linhas, ou seja, pra cada valor na amostra
  eo <- sort(e[i,]) # aqui se trata do residuo ordenado
  e1[i] <- (eo[2]+eo[3])/2 # aqui é o calculo do limite inferior ; pega-se as posições 2 e 3, para se obter o percentil de 2.5
  e2[i] <- (eo[97]+eo[98])/2 } # aqui é o calculo do limite superior ; pega-se as posições 97 e 98, para se obter o percentil de 97.5

##
med <- apply(e,1,mean) # linha de referência (setimo passo da criação do grafico envelope)
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