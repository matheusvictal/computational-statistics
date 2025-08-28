

#######################################################

rm(list=ls(all=TRUE))

library(ggplot2)

#1-) Amostrador de Gibbs



y <- c(98.77, 97.12, 97.34, 98.66, 99.66, 97.22, 98.04, 97.30, 97.15, 98.80,
       98.53, 97.96, 98.05, 97.14, 97.53, 97.50, 98.32, 97.24, 97.70, 98.13,
       99.97, 99.01, 97.54, 98.66, 97.78, 99.34, 98.74, 97.36, 97.66, 97.70,
       97.98, 97.96, 97.26, 99.32, 98.32, 97.29, 98.41, 96.60, 97.37, 98.14,
       97.47, 98.54, 97.25, 97.54, 98.56, 96.89, 97.94, 98.39, 99.80, 99.00, 
       99.06, 97.17, 98.35, 98.40, 97.64, 99.34, 99.90, 96.36, 97.60, 98.33,
       99.46, 97.22, 96.76, 98.66, 98.91, 96.94, 97.19, 98.28, 99.33, 96.98,
       98.36, 98.40, 99.76, 99.34, 98.90, 99.23, 97.78, 99.69, 97.45, 96.94,
       99.20, 98.59, 99.60, 98.32, 98.32, 98.50, 98.83, 97.04, 98.91, 97.24, 
       99.06, 97.79, 98.68, 97.63, 97.53, 98.60, 99.35, 98.95, 98.06, 99.65,
       96.96, 99.15, 97.94, 97.07, 99.13, 97.92, 97.76, 99.21, 98.91, 98.45, 
       99.25, 98.93, 96.03, 98.52, 99.16, 97.64, 98.21, 99.46, 97.42, 99.31,
       98.00, 97.75, 97.99, 99.15, 96.93, 96.79, 98.44, 99.20, 99.31, 97.73)


s <- sum(y)

mu <- function(x2){
  return((1/2+x2*s)/(130*x2+1/100))
}

sigma2 <- function(x2){
  return(1/(130*x2+1/100))
}

k = 65.001

theta <- function(x1){
  return( 1 / ( 0.001 + sum((y-x1)^2)/2 ) )
}


# X1|X2=x2 ~ N(mu,sigma2)

# X2|X1=x1 ~ Gamma(k,theta)

amostrador <- function(iter, x1.inicial, x2.inicial, semente){
  
  set.seed(semente)
  
  R <- iter
  x1 <- x2 <-c()
  
  x1[1] <- x1.inicial
  x2[1] <- x2.inicial
  
  for(j in 1:R) {

    x1[j+1] <-rnorm(1, mu(x2[j]), sqrt(sigma2(x2[j])) )
    x2[j+1] <-rgamma(1, shape = k, scale = theta(x1[j+1]) )
  
    }
  return(list(x1,x2))
}




amostras <- amostrador(iter = 1000, x1.inicial = 10000, x2.inicial = 500, 
                       semente = 2112)


ax1 <- unlist(amostras[1], recursive = TRUE, use.names = TRUE)
ax2 <- unlist(amostras[2], recursive = TRUE, use.names = TRUE)



# Gráficos
par(mfrow =c(1, 2))
plot(ax1[-1], xlab = "Índice", ylab = "X1", type = "l")
plot(ax2[-1], xlab = "Índice", ylab = "X2", type = "l")



summary(ax1[-1])
sd(ax1[-1])
summary(ax2[-1])
sd(ax2[-1])

acf(ax1[-1], lag.max = 10, xlab = "Defasagem", main = "",
    ylab = "Autocor. de X1")
acf(ax2[-1], lag.max = 10, xlab = "Defasagem", main = "",
    ylab = "Autocor. de X2")



#######################################################

rm(list=ls(all=TRUE))

library(ggplot2)


#2-) Teste de permutação

tempos1 <- c(63, 82, 81, 68, 57, 59, 66, 75, 82, 73)

tempos2 <- c(64, 56, 72, 63, 83, 74, 59, 82, 65, 82)

# Teste de igualdade de variabilidade:

## Estatística de teste: s2_1/s2_2;
## Aplicar um teste de permutação (não há suposições sobre a distribuição 
## original);
## Muitas permutações são possíveis: abordagem por método de Monte Carlo.

n1 <- length(tempos1)
n2 <- length(tempos2)

n <- n1 + n2

dados <- c(tempos1,tempos2)

choose(n, n1) # permutações possíveis

F_obs <- var(tempos1)/var(tempos2)
F_obs


plot(ecdf(tempos1), main = "", pch = 20, col = "blue",
     xlab = "Tempos observados",
     ylab = "Distribuição empírica",
     xlim =range(dados))
lines(ecdf(tempos2), pch = 20, col = "red")
rug(tempos1, col = "blue", lwd = 2)
rug(tempos2, col = "red", lwd = 2, side = 3)
legend("topleft",c("Tipo 1", "Tipo 2"), lty = 1, col =c("blue", "red"),
       bty = "n")


R <- 9999
F. <-c()

set.seed(2112)

for(j in 1:R) {
  dadosl <- sample(dados)
  xl <- dadosl[1:n1]
  yl <- dadosl[(n1+1):n]
  F.[j] <- var(xl)/var(yl)
}


hist(F., freq = FALSE, main = "", 
     xlab = "Razão entre variâncias amostrais",
     ylab = "Densidade")
lines(density(F.), col = "magenta", lty = 2, lwd = 2)
points(F_obs, 0, pch = 19, col = "red")
points(1/F_obs, 0, pch = 19, col = "red")
box()




plot(table(F.)/R, col = "gray50", 
     xlab = "Razão entre variâncias amostrais",
     ylab = "Frequência relativa")
points(F_obs, 0, pch = 19, col = "red")
points(1/F_obs, 0, pch = 19, col = "red")
box()


kd <- sum(F. <= F_obs)
ke <- sum(F. >= 1/F_obs)

p1 <- (ke+1)/(R+1)
p2 <- (kd+1)/(R+1)

cat("\n valor-p esquerda aprox. =", (ke+1)/(R+1), "\n")
cat("\n valor-p direita aprox. =", (kd+1)/(R+1), "\n")

p <- 2*min(p1,p2)

cat("\n valor-p aprox. =", p, "\n")


