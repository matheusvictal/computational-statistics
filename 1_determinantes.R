# Aqui utiliza-se de uma solucao computacional elegante para obtencao da 
# distribuicao de frequencias de todos os determinantes possiveis em matrizes
# 2x2 cujos elementos variam de 0 a 9.

prod <- outer(0:9,0:9)
frdet <- table(outer(prod,prod, FUN = '-'))

plot(frdet, xlab = 'Determinante', ylab = 'Frequencia', col = 'firebrick')

# E interessante evitar o uso de lacos (for, while) para que o processamento seja
# mais otimizado