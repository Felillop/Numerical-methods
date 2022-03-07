## PEC 2

install.packages("pracma")
library(pracma)

#Definicion de la matriz L
L = matrix (c(2, 35, -36, -180,
              1,  0,   0,    0,
              0,  1,   0,    0,
              0,  0,   1,    0) , nrow =4, byrow = TRUE )

# Orden de la matriz L
N = dim(L)[1]

# Alternativa en el calculo de los coeficientes del polinomio caracteristico
#Matriz identidad de dimension igual que L
I<-diag(N) 
B1 = L
B1

c1 = -sum(diag(B1))
c1

B2 = L%*%(B1 + c1 * I)
B2

c2 = -(1/2)*sum(diag(B2))
c2

B3 = L%*%(B2+c2 * I)
B3

c3 = -(1/3)*sum(diag(B3))
c3

B4 = L%*%(B3+c3 * I)
B4

c4 = -(1/4)*sum(diag(B4))
c4

#Calculo de los autovalores de L con mayor modulo ----------------------------------------------------------------------------------------

# Como se piden los valores de an y bn para 3<=n<=20, debemos calcular s1, s2, ...,s20

#Construccion de la sucesion de sn desde s1 hasta s20
n = 20  #n denota el número de iterantes que se van a calcular

Ln = diag(N)
sn = 0
for(i in 1:n){
  Ln = Ln%*%L  # Ln multiplicada por L
  sn[i] = sum(diag(Ln)) #Traza de la matriz Ln
}
print(sn)

#calculo de sigma y delta

sigma = c(0, 0) #inicializamos a cero el vector
sigma

for(i in 2:length(sn)){
  sigma[i] <- sn[i]/sn[i-1]
}
print(sigma)

#Inicializamos el vector de esta forma para evitar posibles NA en las posiciones 1 y 2
delta = c(0, 0, 0)
delta
for(i in 3:length(sigma)){
  delta[i] = (sigma[i]-sigma[i-1])/(sigma[i-1]-sigma[i-2])
}

# A continuacion calculamos las sucesiones an y bn.
# Para ellos debemos resolver el siguiente sistema de ecuaciones:
#a[i] + b[i]  = sigma[i] + delta[i]*sigma[i-2]
#a[i]*b[i] = sigma[i-1]*sigma[i-2]*delta[i]

#Creamos las sucesiones auxiliares yn = an + bn, y zn = an*bn

yn = c(0, 0, 0)
zn = c(0, 0, 0)

for (i in 3:n){
  yn[i] = sigma[i] + (delta[i]*sigma[i-2])
}
print(yn)

for (i in 3:n){
  zn[i] = sigma[i-1]*sigma[i-2]*delta[1]
}
print(zn)

#Los valores de an y bn se obtienen al resolver la ecuacion:
#    x^2 - yn[i]x + zn[i]=0 para cada uno de los valores de i. 
#La funcion polyroot devuelve valores complejos y nosotros nos quedamos con la parte real

for (i in 3:n){
  lambda = Re(polyroot(c(zn[i], -yn[i], 1)))
  print(lambda)
}
lambda4 = lambda[2]
lambda3 = lambda[1]

#Calculo de los dos autovalores lambda_1 y lambda_2 de la matriz L ------------------------------------------------------
num = c(1, c1, c2, c3, c4)
den = c(1, -(lambda[1] + lambda[2]), lambda[1]*lambda[2]) #(x-lambda3)*(x-lambda4)

div = deconv(num, den)

print(div[[1]])

#Calculamos las raices del polinomio obtenido
lambda = Re(polyroot(c(rev(div[[1]]))))
lambda2 = lambda[2]
lambda1 = lambda[1]

print(c(lambda1, lambda2, lambda3, lambda4))

