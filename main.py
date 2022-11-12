import math
import matplotlib.pyplot as plt
import numpy as np


dt = 0.01
t = np.arange (0, 0.5+dt, dt)
dx= 0.1
x = np.arange(0,1+dx,dx)
n = x.size
m = t.size
M = np.zeros((m,n))

for i in range (0,m):
    M[i][0] = 2*t[i]
for i in range (0,m):
    M[i][n-1] = -1
for i in range (1,n):
    M[0][i] = x[i]*math.cos(math.pi*x[i])


for i in range (1,n-1):
    f = x[i]*(2-x[i])
    M[1][i]= M[0][i]+f*dt+ (dt*dt)/(2*dx*dx)*(M[0][i-1]-2*M[0][i]+M[0][i+1])

#явная схема
def explicit_MKR(A):
    for k in range(2, m):
        for i in range(1, n - 1):
            A[k][i] = (A[k - 1][i + 1] - 2 * A[k - 1][i] + A[k - 1][i - 1]) *dt *dt / (dx * dx) + 2*A[k - 1][i]-A[k-2][i]
    return A

M_ex = np.copy(M)
M_ex = explicit_MKR(M_ex)

print(' ')
for i in range (0,m):
  print(*M_ex[i])


#неявная схема
M_im = np.copy(M)
print('   ')

A = 1/(dx*dx)
B = (dx*dx+2*dt*dt)/(dx*dx*dt*dt)
C = 1/(dx*dx)

for i in range(2, m):
    F = np.zeros((n,))
    P = np.zeros((n,))
    Q = np.zeros((n,))
    F[0] = 2*M_im[i-1][0]/(dt*dt) - M_im[i-2][0]/(dt*dt)
    P[0] = C/B
    Q[0] = (F[0] + A*M_im[i][0])/B
    for j in range(1, n-1):
        F[j] = 2*M_im[i-1][j]/(dt*dt) - M_im[i-2][j]/(dt*dt)
        P[j] = C / (B - A*P[j-1])
        Q[j] = (F[j] + A*Q[j-1])/(B - A*P[j-1])
    q = n-2
    for k in range(1, n-1):
        M_im[i][q] = P[q]*M_im[i][q+1] + Q[q]
        q = q - 1
for i in range (0,m):
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

plt.plot(x, M_ex[30], label='t1')
plt.plot(x, M_ex[10],label='t2')
plt.plot(x, M_ex[50],label='t3')
plt.title('График при 3-ёх фиксированных t для явной схемы')
plt.xlabel("x")
plt.ylabel("U", rotation='horizontal')
plt.legend()
plt.show()

plt.plot(x, M_im[30], label='t1')
plt.plot(x, M_im[10],label='t2')
plt.plot(x, M_im[50],label='t3')
plt.title('График при 3-ёх фиксированных t для неявной схемы')
plt.xlabel("x")
plt.ylabel("U", rotation='horizontal')
plt.legend()
plt.show()

plt.plot(x, M_ex[25], label='explicit')
plt.plot(x, M_im[25],'k--' , label='implicit')
plt.title('График сравнение явного и неявного метода при фиксированном времени')
plt.xlabel("x")
plt.ylabel("U")
plt.legend()
plt.show()




