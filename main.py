import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

dx = 0.2
x = np.arange(0, 1+dx, dx)
dy = 0.2
y = np.arange(0, 1+dy, dy)
m = x.size
n = y.size
eps = 0.01

U = np.zeros((m,n))

def first(A):
    for i in range (1,m-1):
        for j in range (1, n-1):
            A[i][j]= 1/4 * (A[i-1][j]+A[i+1][j]+A[i][j+1]+A[i][j-1])
    return A
for i in range(0, n):
    U[i][0] = 20*y[n-i-1]
    U[i][n-1] = 30*math.cos(math.pi * y[n-i-1]/2)
    U[0][i] = 30 * math.cos(math.pi * x[i] / 2)
    U[n - 1][i] = 20*x[i]*x[i]



U_k = np.copy(U)
U_kk = first(U)

print(LA.norm(U_kk - U_k))
omega = 0
IT = np.zeros((1,19))
OM = np.zeros((1,19))
for k in range (1,20):
    omega = k*0.1
    iter = 0
    U_k = np.copy(U)
    U_kk = np.zeros((m, n))
    while LA.norm(U_kk - U_k) > eps:
        iter = iter + 1
        U_kk = np.copy(U_k)
        for i in range (1,m-1):
            for j in range (1, n-1):
                U_k[i][j] = U_k[i][j]*(1-omega) + omega/4*(U_k[i-1][j] + U_k[i+1][j] + U_k[i][j-1] + U_k[i][j+1])
    IT[0][k-1] = iter
    OM[0][k-1] = omega



omega = 1.3
U_k = np.copy(U)
U_kk = np.zeros((m, n))
while LA.norm(U_kk - U_k) > eps:
    U_kk = np.copy(U_k)
    for i in range (1,m-1):
        for j in range (1, n-1):
            U_k[i][j] = U_k[i][j]*(1-omega) + omega/4*(U_k[i-1][j] + U_k[i+1][j] + U_k[i][j-1] + U_k[i][j+1])


print (*OM)
print (*IT)
for i in range (0,m):
    print(*U_k[i])

plt.plot(OM, IT, '*')
plt.title("Зависимость итераций от омеги")
plt.show()

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(projection='3d')
yy, xx = np.meshgrid(y,x)
ax.plot_surface(yy, xx, U_k)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('U')
surf = ax.plot_surface(yy,xx, U_k, cmap='inferno')
plt.title('Решение ')
plt.show()







