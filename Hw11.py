import numpy as np
from scipy.integrate import ode

# Problem 3

def rk4(f):
    """Runge–Kutta 4-th order."""
    def solution(t, y, h):
        while True:
            yield y
            k1 = f(t,     y)
            k2 = f(t+h/2, y + h*k1/2)
            k3 = f(t+h/2, y + h*k2/2)
            k4 = f(t+h,   y + h*k3)
            y += h*(k1 + 2*k2 + 2*k3 + k4)/6
            t += h
    return solution

M1 = np.matrix([[-10, 9], [10, -11]])
M2 = np.matrix([[998, 1998], [-999, -1999]])

for M in [M1, M2]:
    print("Using", M)
    for N in [1,10,100,1000,10000]:
        sol = rk4(lambda t, y: M@y)(0, np.array([[1.],[2.]]), 1/N)
        for _ in range(N*10):
            next(sol)
        sol = next(sol)
        print(N, sol)

# Problem 4

Phi = ode(lambda t, y: [y[1], np.sin(np.log(t))/t**2 + y[0]*2/t**2 - y[1]*2/t])
phi = ode(lambda t, y: [y[1], y[0]*2/t**2 - y[1]*2/t])

c2 = (8 - 12*np.sin(np.log(2)) - 4*np.cos(np.log(2)))/70
c1 = 11/10 - c2
v_exact = c1 - 2*c2 - 3/10
v = 1  # guess
for _ in range(100):
    Pv = Phi.set_initial_value([1, v], 1).integrate(2)[0]
    pv = phi.set_initial_value([0, 1], 1).integrate(2)[0]
    print(v - v_exact, (Pv - 2))
    v -= (Pv - 2) / pv
