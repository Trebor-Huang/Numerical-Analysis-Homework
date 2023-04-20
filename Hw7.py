import numpy as np
# Problem 2.

def f(x):
    return np.exp(-x) - np.sin(x)

def steffensen(f, x):
    while (u := f(x+f(x)) - f(x)) != 0:
        x = x - f(x) ** 2 / u
        yield x

print("\nSteffensen, Initial condition x=0:")
solution = steffensen(f, 0)
exact = 0.58853274398186107743245205
for s in solution:
    print("%.5e" % abs(s - exact))

print("\nSteffensen, Initial condition x=3:")
solution = steffensen(f, 3)
exact = 3.09636393241064611562584085
for s in solution:
    print("%.5e" % abs(s - exact))

print("\nIteration for F1:")
exact = 0.58853274398186107743245205
x = 0.5
for _ in range(20):
    x = np.arcsin(np.exp(-x))
    print("%.5e" % abs(x - exact))

print("\nIteration for F2:")
exact = 6.2850492733825866764618695
x = 6
for _ in range(10):
    x = np.exp(-x) - np.sin(x) + x
    print("%.5e" % abs(x - exact))

# Problem 3

def f(x):
    return np.exp(-x**2/2) - np.cos(2*x)/4
def df(x):
    return -x*np.exp(-x**2/2) + np.sin(2*x)/2
def ddf(x):
    return (x**2 - 1) * np.exp(-x**2/2) + np.cos(2*x)

print("\nNewton's method for maxima:")
x = 1
for _ in range(20):
    x = x - df(x) / ddf(x)
    print("%.5e" % x)

print("\nCorrected Newton's method for maxima:")
x = 1
for _ in range(5):
    x = x - 3 * df(x) / ddf(x)
    print("%.5e" % x)
