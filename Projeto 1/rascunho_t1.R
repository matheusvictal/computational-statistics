rm(list=ls(all=TRUE))
library(ggplot2)

#####################################################################

#1-)

# Função densidade de probabilidade de Laplace (g(y): função auxiliar)
d_lap <-function(x) {
  return(0.5* exp(-abs(x)))
}

# Gráfico para d_lap
curve(d_lap, -5, 5, xlab = "y", ylab = "g(y)", col = "blue", lwd = 2)

# Mesmo gráfico com ggplot2
ggplot() + 
  xlim(-7,7) + 
  geom_function(fun = d_lap, colour = "deeppink3",size = 0.7) +
  xlab("y") + 
  ylab("g(y)")


# Função de densidade f(x) da v.a. X da qual se quer obter uma amostra
f_x <-function(x) {
  return(exp(-(abs(x)^3)/3))
}

# Gráfico de f(x)
ggplot() + 
  xlim(-7,7) + 
  geom_function(fun = f_x, colour = "deeppink3",size = 0.7) +
  xlab("x") + 
  ylab("f(x)")


# Função f(x)/g(x), da qual queremos obter o valor máximo para atribuir a M
fgx <-function(x) {
  return(f_x(x)/ d_lap(x))
  }


# Gráfico de f(x)/g(x)
ggplot() + 
  xlim(-7,7) + 
  geom_function(fun = fgx, colour = "deeppink3",size = 0.7) +
  xlab("x") + 
  ylab("f(x)/g(x)")


# Obtenção do valor positivo de x que otimiza a função f(x)/g(x)
x_max <- optimize(fgx, interval=c(-7,7), maximum=TRUE)
x_max$maximum

M <- fgx(x_max$maximum)
M

# Função para Mg(x)
Mg_x <-function(x) {
  return(M*d_lap(x))
}

curve(Mg_x, -5, 5, xlab = "y", ylab = "g(y)", col = "blue", lwd = 2)
curve(f_x, add = TRUE)


# Gráfico de f(x) e Mg(x)
ggplot() + 
  xlim(-7,7) + 
  geom_function(aes(colour = "f(x)"), fun = f_x, size = 0.75) +
  geom_function(aes(colour = "Mg(x)"), fun = Mg_x, size = 0.75) +
  xlab("x") +
  ylab("Densidade")


# Algoritmo para obtenção de a.a. de X via método da aceitação-rejeição

n <- 1000 # tamanho da amostra

set.seed(2112) 

i <- n_t <- 0 # Inicialização do contador i para o tamanho da amostra e 
# do contador n_t para o número de tentativas 

a_X <- c() # vetor auxiliar para armazenamento dos valores observados da v.a. X

while(i<n){
  rej <- TRUE
  
  while(rej){
    n_t <- n_t + 1
    
    # Variável auxiliar
    u <-runif(1)
    y <-ifelse(u<=0.5,log(2*u),-log(2*(1-u)))
    
    if(runif(1) <= fgx(y)/M){
      i <- i + 1
      a_X[i] <- y
      rej <- FALSE
    }
  }
}

cat("\n Tamanho da amostra:", i, "\n Número de tentativas:", n_t)

hist(a_X, freq = FALSE, main = "", xlab = "x", ylab = "Densidade", 
     xlim = c(-3,3))


#####################################################################

#2-)

n <- 1000000
aa <- c()

for(i in 1:n){
  x <- rlnorm(1,0,1)
  log_x <- log(x)
  y <- exp(9 + 3*log_x + rnorm(1,0,1))
  aa[i] <- y/x
}

E <- mean(aa)
E
var_est <- var(aa)/n
ep <- sqrt(var_est)
ep

#####################################################################

#3-)
 

rc <- function(aa){
  n <- length(aa)
  result <- ifelse(mean(aa) >= (sqrt(2/n)*qnorm(0.95)+2),1,0)
  return(result)
}

n_tam <- seq(5,100,by=5) # tamanho de amostra
n <- 50 # simulações para cada tamanho de amostra

taxas_val <- c()

for(i in n_tam){
  taxa <- 0
  for(j in 1:n){
    aa <- rpois(n = i,lambda = 2) 
    taxa <- taxa + rc(aa)
  }
  taxa <- taxa/length(1:n)
  taxas_val[i] <- taxa
}

plot(taxas_val, ylim = c(0,1))
abline(h=0.05, col = "red")



n_sim <- seq(50,3000,by=50)
n <- 10000
taxas_val <- c()


for(i in n_sim){
  taxa <- 0
  for(j in 1:i){
    aa <- rpois(n,lambda = 2) 
    taxa <- taxa + rc(aa)
  }
  taxa <- taxa/i
  taxas_val[i] <- taxa
}
plot(taxas_val, ylim = c(0,0.1))
abline(h=0.05, col = "red")



##################################################################

rm(list=ls(all=TRUE))


rc <- function(aa){
  n <- length(aa)
  rej <- ifelse( mean(aa) >= (sqrt(2/n)*qnorm(0.95)+2), 1, 0)
  return(rej)
}


n_aa <- c(10,100,300,500)
n_sim <- 10000

set.seed(2112)

for(i in n_aa){
  taxa_I <- 0
  for(j in 1:n_sim){
    aa <- rpois(i, lambda = 2)
    taxa_I <- taxa_I + rc(aa)
  }
  taxa_I <- round(taxa_I/n_sim,10)
  cat("\nTaxa de erro tipo I, aa de tamanho",i,":",taxa_I)
}
  
  





rm(list=ls(all=TRUE))


rc <- function(aa){
  n <- length(aa)
  rej <- ifelse( mean(aa) >= (sqrt(2/n)*qnorm(0.95)+2), 0, 1)
  return(rej)
}


n_aa <- c(10,100,300,500)
n_sim <- 1000

lambdas <- seq(2.05,4, by = 0.05)


df <- data.frame()

set.seed(2112)

for(i in n_aa){
  for(j in lambdas){
    taxa_II <- 0 
    for(k in 1:n_sim){
      aa <- rpois(i, lambda = j)
      taxa_II <- taxa_II + rc(aa)
    }
    taxa_II <- taxa_II/n_sim
    poder <- (1-taxa_II)
    
    obs <- c(i,j,poder)
  
    df <- rbind(df,obs)
  }
}


names(df) <- c("tam","lambda", "poder")

plot(x = df[df$tam==10,]$lambda, y = df[df$tam==10,]$poder, 
     main = "Poder para amostra de tamanho 10 (Monte Carlo, 1000 simulações)",
     xlab = "Valores para lambda",
     ylab = "Poder observado")

plot(x = df[df$tam==100,]$lambda, y = df[df$tam==100,]$poder, 
     main = "Poder para amostra de tamanho 100 (Monte Carlo, 1000 simulações)",
     xlab = "Valores para lambda",
     ylab = "Poder observado")

plot(x = df[df$tam==300,]$lambda, y = df[df$tam==300,]$poder, 
     main = "Poder para amostra de tamanho 300 (Monte Carlo, 1000 simulações)",
     xlab = "Valores para lambda",
     ylab = "Poder observado")

plot(x = df[df$tam==500,]$lambda, y = df[df$tam==500,]$poder, 
     main = "Poder para amostra de tamanho 500 (Monte Carlo, 1000 simulações)",
     xlab = "Valores para lambda",
     ylab = "Poder observado")
