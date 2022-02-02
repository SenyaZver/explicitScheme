import scipy.linalg
import numpy as np
from numpy import pi, cos, sin
from pprint import pprint
import matplotlib as mpt
import matplotlib.pyplot as plt
# mpt.use('TkAgg')
# plt.ion()
np.set_printoptions(suppress=True)
def starting_function (x):
    return cos(pi*x) + 2*pi*x

def border_function (t):
    return 1 - sin(2*pi*t)/2 - 3.5*t


def f(t,x):
    return -cos(2*pi*t)*pi - 3.5 + 0.41 * (-sin(pi * x)*pi + 2*pi)

def exact_function (t,x):
    return cos(pi*x) - sin(2*pi*t)/2 + 2*pi*x - 3.5*t

def findMaxDeviation(u, u_exact):
    return np.max(np.abs(u[:][:]-u_exact[:][:]))


def solve (tau, h, a):
    x_grid_size = int(1 / h)
    t_grid_size = int(1 / tau)

    x = np.linspace(0, 1, x_grid_size)
    t = np.linspace(0, 1, t_grid_size)

    u = np.zeros((t_grid_size,x_grid_size))

    u[0, :] = starting_function(x[:])
    u[:, 0] = border_function(t[:])


    b = np.zeros(x_grid_size - 2)

    alpha = (tau*a)/(h)

    for n in range(0,t_grid_size-1):
        for j in range(1,x_grid_size):
            u[n+1][j] = alpha * u[n][j-1] + (1-alpha)*u[n][j] + tau * f(t[n], x[j])


    return u




                                # Calculating u

a = 0.41

h = 0.1
print("tau should be less or equal to", h/a)
print()
tau = float(input())

if (tau > h/a):
    print("wrong tau")
    exit()


u = solve(tau, h, a)


                                # Creating calculated graph
x_grid_size = int(1 / h)
t_grid_size = int(1 / tau)
x = np.linspace(0, 1, x_grid_size)
t = np.linspace(0, 1, t_grid_size)
X, T = np.meshgrid(x, t)


CalculatedGraph = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(T, X, u, 30, cmap='binary')
ax.plot_surface(T, X, u, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('u')
ax.invert_xaxis()
ax.set_title('Calculated Solution')


                               # Creating exact graph
X_exact, T_exact = np.meshgrid(x, t)

U_exact = exact_function(T, X)

ExactGraph = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(T_exact, X_exact, U_exact, 30, cmap='binary')
ax.plot_surface(T_exact, X_exact, U_exact, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('u')
ax.invert_xaxis()
ax.set_title('Exact Solution')

                                # Creating both graphs
BothGraph = plt.figure()
ax_both = plt.axes(projection='3d')
ax_both.plot_surface(T_exact, X_exact, U_exact, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax_both.plot_wireframe(T, X, u, cmap='binary')
ax_both.set_title('Both solutions')
ax_both.set_xlabel('t')
ax_both.set_ylabel('x')
ax_both.set_zlabel('u')
ax_both.invert_xaxis()



print("actual error", findMaxDeviation(u, U_exact))
print("theoretical error", tau + h)
plt.show()