## PEC 1 ##

#Introducir la matriz L
L=matrix(nrow=4,ncol=4, c(2,35,-36,-180,
                          1, 0,  0,   0,
                          0, 1,  0,   0,
                          0, 0,  1,   0),byrow=TRUE)

#Calculo de los valores s1, s2, s3, s4

#Calculamos la traza de las potencias de la matriz L

L2 = L%*%L
L3 = L2%*%L
L4 = L3%*%L

s1 = sum(diag(L))
print(s1)

s2 = sum(diag(L2))
print(s2)

s3 = sum(diag(L3))
print(s3)

s4 = sum(diag(L4))
print(s4)

#Calculamos los valores propios de la matriz L
# o autovalores

lambda <-eigen(L,only.values=TRUE)
print(lambda)

#Los sumamos
lambda[[1]]

s1_1 = sum(lambda[[1]])
print(s1_1)

s2_1 = sum(lambda[[1]] ^2)
print(s2_1)

s3_1 = sum(lambda[[1]] ^3)
print(s3_1)

s4_1 = sum(lambda[[1]] ^4)
print(s4_1)

#Calculo de lo coeficientes del polinomio carcateristico 

c1 = -s1
print(c1)

c2 = (1/2)*(-s2 - c1*s1)
print(c2)

c3 = (1/3)*(-s3 - c1*s2 - c2*s1)
print(c3)

c4 = (1/4)*(-s4 - c1*s3 - c2*s2 - c3*s1)
print(c4)

#El polinomio obtenido es:
cat("x^4", c1, 'x^3', c2, 'x^2', c3, 'x', c4)

