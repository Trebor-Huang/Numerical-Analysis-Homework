#set text(font: "CMU Serif")
= Homework 10

Recall our notation: we have a linear multistep method $ u_(n+k) + sum_(j=0)^(k-1)a_j u_(n+j) = h sum_(j=0)^k b_j f_(n+j). $

== Problem 1

Suppose a method is not stable, then the zero condition is not satisfied. Therefore either

+ We have a root $lambda$ of the characteristic equation such that $|lambda| > 1$. So applying the method to $f(t,y) = 0$ (which satisfies our condition) gives $ u_(n+k) + sum_(j=0)^(k-1) a_j u_(n+j) = 0. $ We have an exact solution $y(t) = 0$, and a numerical solution $u_n = h lambda^n$. As $h$ tends to zero, the error of the initial data goes to zero, however $u_N = h lambda^((b-a)/h)$ tends to infinity. Therefore the method is not convergent on every $f$.
+ We have a root $|lambda| = 1$ with multiplicity greater than one. Then $u_n = sqrt(h) n lambda^n$ is a solution. As $h$ tends to zero the initial error also goes to zero. But $|u_N| = sqrt(h) (b-a)/h = (b-a)/sqrt(h)$ which tends to infinity.

In both cases, we find that for a particular $f$, there is a set of initial conditions $u_i(h)$ that tends to the exact solution as $h -> 0$, but the numerical solution on the whole interval $[a,b]$ does not tend to the exact solution. This proves that the method does not converge.

== Problem 2

It suffices to let $t_n = u_n = 0$. Let $f(t,u) = A + B t + C u + D t^2 + E t u + F u^2 + O(epsilon^3)$. After some arduous Taylor expansion we see that
$ k_1 &= A h \
k_2 &= A h + (B + A C) alpha h^2 + (D + A E + A^2 F) alpha^2 h^3 + O(h^4)\
k_3 &= A h + (B + A C) (1 - alpha) h^2 + (C(B + A C)alpha + D(...) + E(...) + F(...))(1-alpha) h^3 + O(h^4) $
So $u_(n+1) = 1/2(k_1 + k_3) = A h + (1-alpha)/2 (B + A C)h^2 + O(h^3)$. In reality the solution should be (as computed in previous homeworks) $ y_(n+1) = A h + (B + A C)/2 h^2 + O(h^3). $ So the order of truncation error is $1$ for $alpha != 0$, and for $alpha = 0$ after computing one more term we have the order of truncation error is $2$.

For stability, after a perturbation of magnitude less than $Delta$, and suppose $u_n$ is disturbed by $c Delta$, we get
$ |k_1' - k_1| &= h |f(t_n, u_n + delta) + delta_1 - f(t_n u_n)| <= (L c + 1) Delta h \
|k_2' - k_2| &<= (L c + 1) Delta h + L(L c + 1)Delta h^2 \
|k_3' - k_3| &<= (L c + 1) Delta h + L (L c + 1) Delta h^2 + L^2(L c+1) Delta h^3. $
Here $L$ is the Lipschitz constant. Note that being Lipschitz in one variable and uniform in another implies that you can use the same Lipschitz constant once you fix a compact set, which is currently the interval $[a,b]$ on which we wish to solve the ODE. So $u_(n+1)$ is disturbed by no more than $c Delta + (L c + 2) Delta h$ for $h < L/2$. The iteration $c |-> (L h + 1)c + 2 h$ has solution $ c_n = (L h + 1)^n (c_0 + 2/L) - 2/L. $ For fixed interval $[a,b]$ this will be iterated $n = (b-a)/h$ times, giving a bound of the maximum error $ (L h + 1)^((b-a)/h) (c_0 + 2 / L) <= e^(L(b-a))(c_0 + 2/L) $ which is indeed independent of $h$. So the method is stable.


== Problem 3

By Dalhquist's equivalence we only need to find a consistent and unstable method. Consider $ u_(n+2) - 3u_(n+1) + 2 u_n = - h f_n. $ This method is consistent because
$ u_(n+2) - 3u_(n+1) + 2 u_n &= (u_(n+2) - u_(n+1)) - 2(u_(n+1) - u_n) \
&= h (f_n + O(h)) - 2h (f_n + O(h))\
&= -h f_n + O(h^2). $
So the method has order 1. However this method is unstable since the characteristic polynomial $x^2 - 3 x + 2 = 0$ has two roots, $1,2$, which don't satisfy the root condition. Numerically we can test it on $u' = u$. The coding is trivial, see appendix.
#align(center, image("Hw10-Fig3.png", width: 70%))
The numerical solutions diverge very quickly.

== Problem 4

We fix the number of nodes $N$, let $x_k = (2 pi k)/N$, $Delta = (2 pi)/N$ is the spatial step size. We can approximate the derivative $(diff u)/(diff x)$ by the spectral differentiation method, i.e. interpolating $u_k(t) = u(t, x_k)$ with sine waves, and calculating the derivative. However Fourier transforming already leads to the exact solution, so let's just use the approximation $(u_(k+1)(t) - u_(k-1)(t))/(2Delta)$. This turns the PDE into an ODE with $N$ variables.

To solve this ODE, we use the forward Euler method and the Runge--Kutta order 4 method. We compare the error in $ell_oo$ norm. For long term evaluations the numerical solutions tends to explode for large $N$, even though the integrator performs relatively well on single variable ODEs. Convergence is also slow.

== Appendix: Code
Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

#raw(read("Hw10.py"), block: true, lang:"python")
