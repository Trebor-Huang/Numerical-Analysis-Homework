import numpy as np
# Only used to provide exact values of mathematical functions

def estimate(Q):
    """
    Estimate the error by 2[I(h) - I(h/2)].
    This works assuming I(h) ~= I + c h.
    """
    def _Q(f, a, b):
        Ih = Q(f, a, b)
        Ih2 = Q(f, a, (a+b)/2) + Q(f, (a+b)/2, b)
        return Ih2, abs(Ih - Ih2)
    return _Q

@estimate
def simpson(f, a, b):
    return (f(a) + 4*f((a+b)/2) + f(b)) / 6 * (b-a)

_p = 1/np.sqrt(3)
@estimate
def newton(f, a, b):
    p1 = (a * (1+_p) + b * (1-_p))/2
    p2 = (a * (1-_p) + b * (1+_p))/2
    return (f(p1) + f(p2)) / 2 * (b-a)

def Quadrature(Q, f, a, b, eps):
    """
    Adaptive quadrature by bisecting the intervals.
    """
    res, err = Q(f,a,b)
    if err < eps:
        return res
    else:
        return \
            Quadrature(Q, f, a, (a+b)/2, eps/2) \
          + Quadrature(Q, f, (a+b)/2, b, eps/2)

def f1(x):
    return np.log(x) * np.sqrt(x) + np.sin(20*x)
# exact value : np.sin(10)**2/10 - 4/9
f1.left = 1e-8
f1.right = 1
f1.res_cutoff = -0.41484854752239014547381181606161963

def f2(x):
    return np.sin(x*x)
f2.left = 0
f2.right = 10
f2.res_cutoff = 0.583670899929623342157572409285574981

for f in [f1, f2]:
    for eps in [1e-6, 1e-9, 1e-12]:
        print("eps=%.0e\tsimpson: %+.3e\tnewton: %+.3e" %
            (eps,
            Quadrature(simpson, f, f.left, f.right, eps) - f.res_cutoff,
            Quadrature( newton, f, f.left, f.right, eps) - f.res_cutoff))
