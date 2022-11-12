import math
import matplotlib.pyplot as plt
import numpy as np


n = 11
M = np.zeros((n, n))
dt = 0.001
dx = 0.06
t = np.arange(0, 0.011, dt)
x = np.arange(0, 0.62, dx)

for i in range(0, n):
    M[i][0] = 0.435 - 2*t[i]
for i in range(0, n):
    M[i][n-1] = 0.8674
for i in range(1, n-1):
    M[0][i] = math.sin(x[i]+0.45)


#явная схема
def explicit_MKR(A):
    for j in range(1, n):
        for i in range(1, n - 1):
            A[j][i] = (A[j - 1][i + 1] - 2 * A[j - 1][i] + A[j - 1][i - 1]) * dt / (dx * dx) + A[j - 1][i]
    return A


M_ex = np.copy(M)
M_ex = explicit_MKR(M_ex)

print(' ')
for i in range (0,n):
  print(*M_ex[i])


#неявная схема
M_im = np.copy(M)
print('   ')

A = -1/(dx*dx)
B = (dx*dx+2*dt*dt)/(dx*dx*dt*dt)
C = -1/(dx*dx)


for i in range(1, n):
    F = np.zeros((n,))
    P = np.zeros((n,))
    Q = np.zeros((n,))
    F[0] = M_im[i-1][0]
    P[0] = C/B
    Q[0] = (F[0] + A*M_im[i][0])/B
    for j in range(1, n-1):
        F[j] = M_im[i-1][j] / (dt*dt)
        P[j] = C / (B - A*P[j-1])
        Q[j] = (F[j] + A * Q[j-1]) / (B - A*P[j-1])
    q = n-2
    for k in range(1, n-1):
        M_im[i][q] = P[q]*M_im[i][q+1] + Q[q]
        q = q - 1

for i in range (0,n):
    print(*M_im[i])

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(projection='3d')
tt, xx = np.meshgrid(x, t)
ax.plot_surface(tt, xx, M_ex)
ax.set_xlabel('x, с')
ax.set_ylabel('t, м')
ax.set_zlabel('U')
surf = ax.plot_surface(tt, xx, M_ex, cmap='inferno')
plt.title('Явная схема решения в зависимости от x,t')
plt.show()


fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(projection='3d')
tt, xx = np.meshgrid(x, t)
ax.plot_surface(tt, xx, M_im)
ax.set_xlabel('x, с')
ax.set_ylabel('t, м')
ax.set_zlabel('U')
surf = ax.plot_surface(tt, xx, M_im, cmap='inferno')
plt.title('Неявная схема решения в зависимости от x,t')
plt.show()

plt.plot(x, M_ex[5], label='t1')
plt.plot(x, M_ex[1],label='t2')
plt.plot(x, M_ex[3],label='t3')
plt.title('График при 3-ёх фиксированных t для явной схемы')
plt.xlabel("x")
plt.ylabel("U", rotation='horizontal')
plt.legend()
plt.show()

plt.plot(x, M_im[3], label='t1')
plt.plot(x, M_im[1],label='t2')
plt.plot(x, M_im[5],label='t3')
plt.title('График при 3-ёх фиксированных t для неявной схемы')
plt.xlabel("x")
plt.ylabel("U", rotation='horizontal')
plt.legend()
plt.show()

plt.plot(x, M_ex[5], label='explicit')
plt.plot(x, M_im[5],'k--' , label='implicit')
plt.title('График сравнение явного и неявного метода при фиксированном времени')
plt.xlabel("x")
plt.ylabel("U")
plt.legend()
plt.show()


