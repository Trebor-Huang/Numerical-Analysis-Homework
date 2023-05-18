#set text(font: "CMU Serif")
= Homework 8

== Problem 1

For the Jacobi method
$ x_(n+1) = D^(-1) b - D^(-1) (L + U) x_n, $
we need to confirm that
$ || D^(-1)(L+U) || = limsup_x (|D^(-1)(L+U)x|) / (|x|) $
is less than one. Here we use the taxicab norm $|x| = sum_i |x_i|$. This is true, because we can calculate the norm of the resulting vector
$ sum_i |y_i| = sum_i |1/(a_(i i)) sum_(i != j) a_(i j) x_j| <= sum_(i != j) |a_(i j)|/|a_(i i)| |x_j| <= M sum_j |x_j|, $ where $M = max_(i) sum_(j != i) |a_(i j)| / |a_(i i)| < 1$. This means it is a contracting map, and thus has a unique fixpoint.

For the Gauss--Seidel method $ x_(n+1) = (D+L)^(-1)b - (D+L)^(-1)U x_n , $ we need a lemma first.

/ Lemma.: If the spectral radius of a matrix $A$ is less than $1$, then it is a contracting map in the taxicab norm when the basis is chosen suitably.
_Proof_. We prove the complex case first. Consider the Jordan normal form of $A$. We have a set of eigenvectors $v_k$ such that $A v_k = lambda_k v_k$ or $A v_k = lambda_k v_k + v_(k+1)$. We scale $u_k = epsilon^(-k) v_k$. This turns the latter case into $A u_k = lambda_k u_k + epsilon u_(k+1)$ for arbitrary $epsilon$. We choose this $epsilon$ so that $max_k|lambda_k| + epsilon < 1$. Now we consider the taxicab norm in the basis ${u_k}$, i.e. the vector $x = sum_i c_i v_i$ has norm $sum_i |c_i|$. Now for $A x = sum_i d_i v_i$, either $d_i = lambda_i c_i$, or $d_i = lambda_i c_i + epsilon c_(i-1)$. So the norm $ sum_i |d_i| <= sum_i |lambda_i| |c_i| + epsilon |c_i| <= (max_k|lambda_k| + epsilon) sum_i |c_i|. $ This proves the claim for complex cases. For the real case we take the real part of everything, and the results transfer through.

Now we need only to estimate the spectral radius of $(D+L)^(-1) U$. So we consider the largest solution in absolute value of $det((D+L)^(-1) U - lambda I) = 0$, or equivalently $ det(U - lambda(D + L)) = 0. $ If this is true, then there is a vector $v$ such that $U v = lambda (D + L) v$. Rearranging we have $ lambda a_(i i) v_i = sum_(j = i+1)^n a_(i j) v_j - lambda sum_(j=1)^(i-1) a_(i j) v_j. $ Therefore $ |lambda| |a_(i i)| |v_i| <= sum_(j = i+1)^n |a_(i j)| |v_j| + |lambda| sum_(j=1)^(i-1) |a_(i j)| |v_j|. $ Consider the largest component of $v$, i.e. index $i$ such that $|v_i| >= |v_j|$ for all other $j$. This gives $ |lambda| |a_(i i)| <= sum_(j = i+1)^n |a_(i j)| + |lambda| sum_(j=1)^(i-1) |a_(i j)| $ or equivalently $ |lambda| <= (sum_(j = i+1)^n |a_(i j)|) / (|a_(i i)| - sum_(j=1)^(i-1) |a_(i j)|). $ Since the matrix is diagonally dominant, we have $|a_(i i)| > sum_(j=1)^(i-1) |a_(i j)| + sum_(j = i+1)^n |a_(i j)|$, therefore $|lambda| < 1$. This proves the result.

== Problem 2

For (i), the "only if" part is definitely false, because if we happen to choose the actual solution, then the iteration will converge no matter what. For the "if" part, we need to see $(I-theta A)$ is a contracting map. By the previous lemma, we see that we need the spectral radius of $(I - theta A)$ to be less than $1$. Suppose $lambda$ is an eigenvalue, i.e. $det(I - theta A - lambda I) = 0$. Then $det(A - (1-lambda) / theta I) = 0$ too, implying $(1-lambda)/theta$ is an eigenvalue of $A$. Therefore $0 < (1-lambda)/theta <= r_sigma(A)$ by symmetry and positive definacy. So $1 - theta r_sigma(A) <= lambda < 1$. As long as $0 < theta r_sigma(A) < 2$, we have our result.

For (ii), consider the largest $rho_"max"$ and smallest eigenvalue $rho_"min"$ of $A$, we have $ 1 - theta rho_"max" <= lambda <= 1 - theta rho_"min", $ therefore we should choose $theta$ such that $theta rho_"max" - 1 = 1 - theta rho_"min"$, which means $ theta = 2/(rho_"max" - rho_"min"). $

== Problem 3

Consider the map $U$ defined as
$ u |-> integral_a^- k(-, s, u(s)) dif s + f. $
This is a map $C[a,b] -> C[a,b]$. Consider the uniform norm $||u|| = sup_x |u(x)|$. Under this norm, $ sup_t | integral_a^t k(t, s, u_1(s)) - k(t,s,u_2(s)) dif s | &<= integral_a^b |k(t,s,u_1(s)) - k(t,s,u_2(s))| dif s \ &<= M integral_a^b |u_1(s) - u_2(s)| dif s \ &<= M (b-a) dot.c sup_t |u_1(s) - u_2(s)|. $ Therefore $||U u|| <= M(b-a) ||u||$. For a small enough interval, this is a contracting map, and therefore has a unique fixpoint. Now note that a Volterra equation can be split up on subintervals, each using a different $f$ but the same $k$: The first interval uses the original $f$, and suppose we have solved the equation on the interval $[a, c]$, then on the next interval $[c,d]$ we have $ u(t) = integral_c^t k(t,s,u(s)) dif s + f(t) + underbrace(integral_a^c k(t,s,u(s)), "known constant"). $ Therefore we just use $f + "the constant"$ as our new $f$. We can piece the solutions on each subinterval to get the desired unique solution.

== Problem 4

The error for $y_(n+1)$ is $ R_(n+1) = |integral_(t_(n-1))^(t_(n+1)) f(s, y(s)) dif s  - 2h f(t_n, y(t_(n)))| <= integral_(t_(n-1))^(t_(n+1)) |f(s, y(s)) - f(t_n, y(t_n))| dif s. $ By Lipschitz continuity, the integrand is less than $h^2 (K + L M)$. This results in an error twice as large as the forward Euler method.

We try out the algorithm on three different ODEs. The first one is $ y' = -y $ with exact solution $y = c e^(-t)$. For $h = 1$ or $h = 0.1$, the numerical solution quickly diverges. For $h = 0.01$, the solution oscillates around $ plus.minus 10^(-3)$. But for $h=10^(-3)$ the error again raises to around $10^(-1)$ for the first one thousand steps. The oscillation is because, when $y_(n-1)$ underestimates and $y_n$ overestimates, both errors cause a underestimation of $y_(n+1)$. Therefore the next step will cause an overestimation, _ad infinitum_.

The second one is $ y' = i y $ with solution $y = exp(i t).$ However the numerical solution has error around $1$ no matter how small the step size is. The solution effectively rotates at a different rate. This is because the secant estimation of the derivative causes an error of the factor $sin(h)/h$.

The final one is $ y' = t - y + 1 $ with a solution $y = e^(-t)+t$. The error _increases_ with the step size at a quadratic rate. The solution largely relies on cancellation between $t$ and $y$.

== Appendix: Code
Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

#raw(read("Hw8.py"), block: true, lang:"python")
