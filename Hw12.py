import numpy as np
import matplotlib.pyplot as plt
eps = 1e-12

def gaussian(mat, vec=None):
    """Partial pivoting Gaussian elimination.
    Will overwrite the input."""
    if mat.size == 0: return
    # Pivot
    pivot = np.argmax(np.abs(mat[:,0]))
    if pivot != 0:
        mat[[0,pivot],:] = mat[[pivot,0],:]
        if vec is not None:
            vec[[0,pivot]] = vec[[pivot,0]]
    # Eliminate
    if np.abs(mat[0,0]) < eps:
        # It's already done, move on
        gaussian(mat[:,1:], vec)
        return
    if vec is not None:
        vec[1:] -= mat[1:,0] * vec[0] / mat[0,0]
    mat[1:,:] -= np.outer(mat[1:,0]/mat[0,0], mat[0,:])
    gaussian(mat[1:,1:], None if vec is None else vec[1:])

def solve(umat, vec):
    """Finish up the solving procedure."""
    if umat.size == 0:
        return
    if np.abs(umat[-1,-1]) < eps:
        raise ZeroDivisionError
    vec[-1] /= umat[-1, -1]
    umat[-1,-1] = 1
    vec[:-1] -= vec[-1] * umat[:-1,-1]
    solve(umat[:-1,:-1],vec[:-1])


def hilbert(n):
    return np.array([[1./(i+j+1) for i in range(n)] for j in range(n)])

def bad(n):
    return np.array([[1. if i == j or i == n-1 else 0. if i > j else -1. for i in range(n)] for j in range(n)])

for n in range(5,21):
    print("%.3e" % np.linalg.cond(hilbert(n), p=float("inf")))

print()

list_error = []
for n in range(5,51):
    errors = 0
    H = bad(n)
    for _ in range(50):
        x = np.random.randn(n)
        y = H @ x
        gaussian(H, y)
        solve(H, y)
        errors += np.max(np.abs(x - y))
    list_error.append(errors/50)

plt.semilogy(range(5,51), list_error)
plt.show()
