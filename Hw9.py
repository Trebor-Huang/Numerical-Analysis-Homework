import numpy as np  # to provide the accurate answer to compare against
## Problem 2
# u1 = u + 1/6 h[4f(t,u) + 2f(t1,u1) + hf(t,u)]

def sol(h, u):
    # f(t,y) = t + y
    t = 0
    while True:
        yield u
        u = (u + ((6+h)*t + (4+h)*u + 2*h) * h / 6)/(1-h/3)
        t += h

for n in [1, 10, 100, 1000, 10000, 100000, 1000000]:
    s = sol(1/n, 0)
    for _ in range(n): next(s)
    q = next(s)
    print("%.3e" % abs(q - np.e + 2))

print("====")

def fix(h, u):
    c = u - (4+h)*h/6 * np.cos(u)
    u, u1 = c - h/3 * np.cos(u), u
    while abs(u - u1) > 1e-15:
        u, u1 = c - h/3 * np.cos(u), u
    return u

def sol(h, u):
    # f(t,y) = -cos(y)
    while True:
        yield u
        u = fix(h, u)

for n in [1, 10, 100, 1000, 10000, 100000, 1000000]:
    s = sol(1/n, 2*np.arctan(np.tanh(1)))
    for _ in range(n): next(s)
    q = next(s)
    print("%.3e" % abs(q - 2*np.arctan(np.tanh(1/2))))
