import numpy as np
import scipy.special as sp
import scipy.linalg as la
import matplotlib.pyplot as plt

problem = ""

if problem == "2.1":
    x = np.linspace(1,2)
    y1 = np.exp(x)
    y2 = np.sin(x)
    y3 = sp.gamma(x)
    y = 1/x

    c, res, rnk, sing = np.linalg.lstsq(np.stack([y1,y2,y3],axis=-1), y, rcond=None)
    print(c)

    y0 = c[0] * y1 + c[1] * y2 + c[2] * y3
    plt.plot(x, y, label="Accurate")
    plt.plot(x, y0, label="Fit")
    plt.legend()
    plt.show()

if problem == "2.2":
    x = np.linspace(0,1)
    y1 = np.exp(x)
    y2 = np.sin(x)
    y3 = sp.gamma(x)
    y = 1/x

    dy = y - y3
    dy[0] = np.euler_gamma

    c, res, rnk, sing = np.linalg.lstsq(np.stack([y1,y2],axis=-1), dy, rcond=None)
    print(c)

    y0 = c[0] * y1 + c[1] * y2 + 1 * y3
    plt.plot(x, y, label="Accurate")
    plt.plot(x, y0, label="Fit")
    plt.legend()
    plt.show()

# In fact all the built-ins used are compliant with MatLab.

t = np.linspace(0,1)
A = np.vander(t, 12)
b = np.cos(4*t)

# standard answer

x_standard = np.polyfit(t, b, 11)

# i

R = la.cholesky(A.T @ A) # RT * R = AT * A
w = la.solve_triangular(R.T, A.T @ b, lower=True)
x_i = la.solve_triangular(R, w, lower=False)
print(np.max(np.abs((x_i - x_standard)))/np.max(np.abs(x_standard)))

# ii

def qr_mgs(A):
    V = np.copy(A)
    m,n = np.shape(V)
    R = np.zeros((n,n))
    Q = np.zeros((m,n))
    for i in range(n):
        R[i,i] = la.norm(V[:,i])
        Q[:,i] = V[:,i] / R[i,i]
        for j in range(i,n):
            R[i,j] = np.dot(Q[:,i].conj(), V[:,j])
            V[:,j] -= R[i,j] * Q[:,i]
    return Q, R

Q, R = qr_mgs(A)
# print("QR error: ", np.max(np.abs(Q@R - A)))
x_ii = la.solve_triangular(R, Q.T.conj() @ b)
print(np.max(np.abs(x_ii - x_standard))/np.max(np.abs(x_standard)))

# iii

def qr_ht(A, b):
    R = np.copy(A)
    m,n = np.shape(A)
    Qb = np.copy(b)
    for k in range(n):
        v = np.copy(R[k:, k])
        v[0] += (1 if v[0] >= 0 else -1) * la.norm(v)
        v /= la.norm(v)
        R[k:, k:] -= 2 * np.outer(v, (v.T @ R[k:, k:]))
        # R[k+1:, k] = 0
        # actually these are never used, so we can just leave them
        Qb[k:] -= 2 * v * np.dot(v, Qb[k:])
    return Qb[:n], R[:n, :n]

Qb, R = qr_ht(A, b)
x_iii = la.solve_triangular(R, Qb)
print(np.max(np.abs((x_iii - x_standard)))/np.max(np.abs(x_standard)))

# iv

Q, R, P = la.qr(A, mode='economic', pivoting=True) # AP = QR
x_iv = la.solve_triangular(R, Q.T.conj() @ b)  # remember to pivot
print(np.max(np.abs((x_iv - x_standard[P])))/np.max(np.abs(x_standard)))

# v

x_v = la.pinv(A) @ b
print(np.max(np.abs((x_v - x_standard)))/np.max(np.abs(x_standard)))

# vi

U, s, Vh = la.svd(A, full_matrices=False)
x_vi = Vh.T.conj() @ ((U.T.conj() @ b) / s)
print(np.max(np.abs((x_vi - x_standard)))/np.max(np.abs(x_standard)))

