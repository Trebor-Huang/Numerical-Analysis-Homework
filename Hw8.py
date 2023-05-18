import numpy as np
def cEuler(h, y0, y1, f):
    """
    Input step size, the first two values,
    and a binary function. Outputs an iterable
    for the solution of y' = f(t, y).
    """
    yield y0
    yield y1
    t = h
    while True:
        y0, y1 = y1, y0 + 2*h*f(t, y1)
        t += h
        yield y1

def test(h, f, sol, n=100):
    num = cEuler(h, sol(0), sol(h), f)
    for i in range(n):
        q = next(num)
        r = sol(i/h)
        print("%+5e" % abs(q-r))


test(0.1, lambda t,y: t-y+1, lambda t: np.exp(-t) + t, n=10)
