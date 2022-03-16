## Solucion de sistema de ecuaciones mediante metodos directos

# Resolucion del sistema mediante factorizacion LU
install.packages("pracma")
library(pracma)

require(pracma)
# Definicion del  sistema Ax = b
A = matrix (c(1, 2,-2,-1, 1,
              0,-3, 1, 2, 2,
              2,-2,-4, 2, 7,
              3,-3,-1, 5, 9,
              0, 0, 0, 0, 2) , nrow =5, byrow = TRUE )

b = c(-2,
       3,
      -3,
      13,
       2)
  
#Factorizacion LU

mlu = lu(A, scheme='ijk') # metodo de Doolittle

L = mlu$L
U = mlu$U

print(L)
print(U)

# Resolucion del sistema a partir de L y U

y = solve(L, b) #Ly = b                                #Los dos sistemas dan el mismo resultado#
x = solve(U, y) #Ux = y

print(x)

#Comprobamos que obtenemos el mismo resultado, resolviendo de forma directa con solve
x_ = solve(A, b) #Ax_ = b

print(x_)

#==================================================APARTADO 2=================================================================================

## Solucion de sistema de ecuaciones mediante metodos iterativos
# Empleamos el comando de R itersolve(A, b, x0, tol, method), cuyos argumentos son:
# A: la matriz del sistema
# b: el lado derecho del sistema
# x0: aproximacio inicial
# tol: condici?n de parada en base al error cometido
# method: metodo a utilizar ("Gauss-Seidel" o "Jacobi")

require(pracma)


# Resolucion del sistema mediante Jacobi

n = 40 #20 #40, #dimension del sistema

A = diag(n)

print(A)

for (i in 1:n-1) { A[i,i+1]<-1/2}
for (i in 2:n) { A[i,i-1]<-1/2}

b = rep(0,n)
b[1]<-1/2

# Primera iteracion

D = diag(diag(A))
L = -tril(A, -1)
U = -triu(A, 1)

J = inv(D) %*% (L+U)
c = inv(D) %*% b

print(J)
print(c)

  
#Calculamos los autovalores de J

lambda <-eigen(J,only.values=TRUE)

print(lambda)

#calculamos el radio espectral (el máximo de los valores absolutos de los autovalores) # https://www.ugr.es/~anpalom/practica7.html#17
max_autovalor = max(abs(J))

print(max_autovalor)

# Primera iteracion

x0 = rep(0,n)
x1 = J%*%x0 + c

print(x1)


# Resolucion iterativa
sol = itersolve(A, b, x0, tol = 1e-6, method = "Jacobi") #solucion con un error de aproximacion maximo
print(sol)
sol1 = itersolve(A, b, x0, nmax = 80, method = "Jacobi") #Solución con un numero maximo de iteraciones
print(sol1)

#calcular el error relativo
x_ = solve (A, b)
dif_J = sol$x -x_
error_J = norm(dif_J, "2")
error_rel_J = error_J/norm(x_, "2")
print(error_rel_J)

# Resolucion del sistema mediante Gauss-Seidel

n =  40 #20 #40 #dimension del sistema
  
A = diag(n)
for(i in 1:n-1) { A[i,i+1]<-1/2}
for(i in 2:n) { A[i,i-1]<-1/2}

b = rep(0,n)
b[1]<-1/2

print(b)

# Primera iteracion
D = diag(diag(A))
L = -tril(A, -1)
M = D - L
U = -triu(A, 1)

G = inv (M)%*%U         #M^{-1}*U
d = inv (M)%*%b         #M^{-1}*b

print (G)
print (d)

x0 = rep(0,n)
x1 = G %*% x0 + d

print(x0)
print(x1)

# Resolucion iterativa
sol = itersolve(A, b, x0,  tol = 1e-6, method = "Gauss-Seidel") #solucion con un error de aproximacion maximo
print(sol)
sol1 = itersolve(A, b, x0, nmax = 80 , method = "Gauss-Seidel") #Solución con un numero maximo de iteraciones
print(sol1)

#calculo del error relativo     Se realiza en el último momento para incluir x_ 


#Comprobamos que obtenemos el mismo resultado, resolviendo de forma directa con solve
x_ = solve(A, b)

dif_GS = sol$x -x_

print(dif_GS)

# Error relativo
error_GS = norm(dif_GS, "2")

error_rel_GS = error_GS/norm (x_, "2")

print(error_rel_GS)

