#set text(font: "CMU Serif")
= Homework 11

Recall the general formula for Runge--Kutta methods (to fix all the symbol names):
$ u_(n+1) = u_n + h sum_(i=1)^s b_i K_i $
$ K_i = f(t_n + c_i h, u_n + h sum_(i=1)^s a_(i j) K_j). $

== Problem 1

$f(t, y) = lambda y$. Therefore
$ arrow(K) = lambda u_n + h lambda A arrow(K) $
The solution is then
$ arrow(K) = lambda u_n (1 - h lambda A)^(-1) arrow(1). $
The entry of the inverse of a matrix is always of the form $p(h lambda)/det(1-h lambda A)$, where $p$ comes from the matrix adjugate, and is a polynomial of degree $(s-1)$. $det(1+h lambda A)$ is also a polynomial in $h lambda$ of degree $s$. Therefore $ h arrow(b) dot arrow(K) = u_n r(h lambda), $ where $r$ is a rational function with degree $s\/s$. Therefore $u_n = c dot [r(h lambda)]^n$, where $c$ is a fixed constant.

When the method is explicit, $A$ is a lower triangular matrix with main diagonal zero, therefore the determinant (which is the denominator) of $(1 - h lambda A)$ is $1$. This makes $r$ a degree $s$ polynomial.

== Problem 2

We have $u_n = c dot.c [r(h lambda)]^n$. With the same initial condition, the exact solution is $y_n = c [exp(h lambda)]^n$. For it to be of order $s$, we need $ r(h lambda) - exp(h lambda) = O(h^(s + 1)). $ Note that $r$ is a degree $s$ polynomial. By Taylor expansion this leaves only the possibility $ r(z) = 1 + z + z^2/(2!) + dots.c + z^s/(s!). $

== Problem 3

For (i), we compute the eigenvalues $(lambda + 10) (lambda + 11) - 90 = 0$, getting $lambda = -2$ or $lambda = -20$. Therefore the stiffness quotient is $20/2 = 10$. For (ii) we similarly compute the eigenvalues to be $lambda = -1$ and $lambda = -1000$, so the stiffness quotient is $1000$.

We need the step size to be less than around $1/20$ and $1/1000$ to prevent divergence, because $h lambda$ needs to fall in the absolute stability region.

We verify this by evaluating the differential equations on $t in [0,10]$. Starting from the initial condition $vec(1,2)$, the first equation overflows with $h = 1$, but resulted in the accurate solution for $h = 10, 100, 1000, 10000$. For the second equation, $h = 1, 10, 100$ all give NaN as an answer, while $h = 1000, 10000$ are approximately correct. This agrees with the analysis.

== Problem 4

With the shooting method, we define $Phi(u, v, t)$ to be the solution $y(t)$ with the initial condition $y(1) = u, y'(1) = v$. Therefore
$ Phi_(t t)(u, v, t) + 2/t Phi_t(u, v, t) - 2/t^2 Phi(u, v, t) = sin(ln t)/t^2. $
So $phi(u, v, t) = Phi_v(u, v, t)$ satisfies the initial value problem
$ phi_(t t)(u, v, t) + 2/t phi_t(u, v, t) - 2/t^2 phi(u, v, t) = 0 $
where $phi(u, v, 1) = 0, phi_t(u, v, 1) = 1$. We have $u = 1$, and we just need to iterate $ v_(k+1) = v_k - (Phi(u, v_k, 2) - 2)/phi(u, v_k, 2). $
Using this method, the errors of $Phi(u, v_k, 2)$ is $4.8 times 10^(-2), 7.6 times 10^(-7), 2.7 times 10^(-15)$ for the first three iterations. After that the error remains at around $10^(-15)$. The final solution is $v approx 0.91762$, which converges with a fixed error related to the integrator tolarance.

== Appendix: Code
Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

#raw(read("Hw11.py"), block: true, lang:"python")
