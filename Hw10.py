import numpy as np
import matplotlib.pyplot as plt

problem = 4  # Enter problem number

if problem == 3:

    def integrate(h):
        u1, u0 = np.exp(h), 1
        while True:
            yield u0
            u1, u0 = 3*u1 - 2*u0 - h * u0, u1

    # The accurate solution
    time = np.linspace(0,1)
    plt.plot(time, np.exp(time), label="Accurate")

    # Increasing N
    for N in [10, 100, 1000]:
        res = integrate(1/N)
        res = np.array([next(res) for _ in range(N)])
        time = np.arange(N)/N
        plt.plot(time, res, label=f"h = 1/{N}")
    plt.axis((0,1,-3,3))
    plt.legend(loc="lower right")
    plt.show()

if problem == 4:
    def forward_euler(f):
        """
        Forward Euler method to solve the equation
          y' = f(t, y).
        Returns an iterator.
        """
        def solution(t, y, h):
            while True:
                yield y
                y += h * f(t, y)
                t += h
        return solution

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

    N = 100
    xs = 2*np.pi*np.arange(N)/N

    def f(t, y):
        return -(np.roll(y, 1) - np.roll(y, -1)) / (4 * np.pi) * N

    def accurate(t):
        return np.sin(xs - t)

    def get_error(E, solution, accurate):
        for _ in range(2**E): next(solution)
        sol = next(solution)
        return np.max(np.abs(sol - accurate))

    test = rk4(f)(0, np.sin(xs), 1/10000)
    for i in range(100):
        r = next(test)
        print(np.max(np.abs(r - accurate(i/10000))))

    errFE = []
    errRK = []
    for E in range(2,10):
        errFE.append(get_error(E-2, forward_euler(f)(0,np.sin(xs),2**(-E)), accurate(0.25)))
        errRK.append(get_error(E-2, rk4(f)(0,np.sin(xs),2**(-E)), accurate(0.25)))

    print(errFE, errRK)
