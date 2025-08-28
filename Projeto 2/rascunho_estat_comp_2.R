rm(list=ls(all=TRUE))

setwd("C:/Users/mathe/OneDrive/Área de Trabalho/1sem 2021/Estatística Computacional/Trabalho 2")

data <- read.csv("dados_IBGE_trat_2018.csv", 
                 header = TRUE,
                 encoding = "UTF-8", dec = ',')


v <- c(1,2,3,4,5,1,2,3,4)

sum(abs(outer(v,v, FUN = "-")))



# Gini coefficient
Gini <- function(v){
  
  mu <- mean(v)
  n <- length(v)
  
  m_sub <- outer(v,v, FUN = '-')
  mod_m_sub <- abs(m_sub)
  somatorio <- sum(mod_m_sub)
  
  gini <- somatorio/(2*n^2*mu)
  
  return(gini)
}

# Population Gini coefficient unbiased estimator 
Gini_est <- function(v){
  
  mu <- mean(v)
  n <- length(v)
  
  m_sub <- outer(v,v, FUN = '-')
  mod_m_sub <- abs(m_sub)
  somatorio <- sum(mod_m_sub)
  
  gini <- somatorio/(2*n^2*mu)
  
  est <- (n/(n-1))*gini
  
  return(est)
}



#Tests

ineq <- c(rep(0,100), 1)
eq <- c(17,17,17,17,17,17,17,17,17,17)

Gini(ineq)
Gini(eq)




# Population Gini Coefficient:

print(Gini(data$PIBpc))



## Bootstrap
set.seed(2112)


n <- length(data$PIBpc)
B <- 8000  
gs <- c()

for (b in 1:B) {
  gs[b] <- Gini(sample(data$PIBpc, n, replace = TRUE))
}

hist(gs, freq = FALSE, main = "", xlab = expression("G"^"*"),
     ylab = "Densidade", col = "lightyellow") 
lines(density(gs),col = "blue", lty = 2, lwd = 2)

box()

# Intervalos de confiança
conf <- 0.95
problu <- c(1 - conf, 1 + conf) / 2
cat("\n Intervalo percentil de", 100 * conf, 
    "%: [",quantile(gs, probs = problu, type = 6), "]")


G_pop <- Gini(data$PIBpc)
G_pop










PIBpc <- data$PIBpc 
n <- length(PIBpc)


## Bootstrap iterado
ga <-  Gini(PIBpc) # Estimativa
set.seed(2112)
B <- 1000  # N de amostras nível 1
B2 <- 50  # N de amostras nível 2 
te <- c()  # tb*

system.time({
  # Nível 1
  ge <- c()
  for (b in 1:B) {
    # amostra bootstrap
    ab <- sample(PIBpc, n, replace = TRUE)
    # Estimativa (*)
    ge[b] <- Gini(ab)
    # Nível 2 
    gee <- c()
    for (j in 1:B2) {
      # Estimativa (**)
      gee[j] <- Gini(sample(ab, n, replace = TRUE))
      }
    te[b] <- (ge[b] - ga) / sd(gee)
  }
})





hist(te, freq = FALSE, main = "", xlab = expression(t^"*"),ylab = "Densidade", 
      col = "lightgreen")
lines(density(te), col = "blue", lty = 2, lwd = 2)
box()




epb <- sd(ge) # erro padrão bootstrap 

# Intervalos de confiança
conf <- 0.95
problu <- c(1 - conf, 1 + conf) / 2
cat("\n Intervalo percentil de", 100 * conf, 
    "%: [",quantile(ge, probs = problu, type = 6), "]")


q12 <- quantile(te, probs = problu, type = 6)
ictb <- ga - q12[2:1]*epb
cat("\n B2:", B2, "\n Intervalo t bootstrap de",
                              100 * conf,"%: [", ictb, "]")


# Visualização
par(mai = c(1, 1, 0.3, 0.2))

hist(ge, freq = FALSE, main = "", xlab = expression(G^"*"),
     ylab = "Densidade", cex.axis = 1.5, cex.lab = 1.5)
lines(density(ge), col = "blue", lwd = 2)
arrows(ictb[1], 0, ictb[2], 0, col = "red", code = 3,
       angle = 90, lwd = 2, length = 0.1)
box()


data <- read.csv2("dadosIBGE_f.csv", 
                 header = TRUE,
                 encoding = "UTF-8", sep = ",") # dados dos municípios de SP, IBGE
