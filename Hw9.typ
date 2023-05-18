#set text(font: "CMU Serif")
= Homework 9

== Problem 1

*Discrete case*. All indices and summations start at zero unless otherwise noted. Given non-negative numbers $g_0, k_s, p_s$ and $phi_n <= g_0 + sum_(s < n) (k_s phi_s + p_s),$
we need to prove $ phi_n <= (g_0 + sum_(s < n) p_s)exp(sum_(s < n)k_s). $ We can let $q_n = g_0 + sum_(s<n) p_s$ be a non-decreasing sequence. This simplifies the inequality to proving $phi_n <= q_n exp(sum_(s<n)k_s)$ given $phi_n <= q_n + sum_(s<n) k_s phi_s$. By induction, we only need to check whether
$ q_n + sum_(s < n) k_s q_s exp(sum_(t < s) k_s) <= q_n exp(sum_(s<n) k_s). $ Since $q_n$ is non-decreasing, it suffices to prove $1 + sum_(s<n) k_s exp(sum_(t<s) k_s) <= exp(sum_(s<n) k_s)$. We proceed by induction again. Suppose this inequality holds for $n$ numbers, we apply it on $k_1, ..., k_n$, giving us $ 1 + sum_(1<=s<=n) k_s exp(sum_(1<=t<s) k_s) <= exp(sum_(1<=s<=n) k_s). $
Multiplying by $exp(k_0)$, we get
$ exp(sum_(s<=n) k_s) &>= exp(k_0) + sum_(1<=s<=n) k_s exp(sum_(t<s) k_s) \ &>= 1 + k_0 + sum_(1<=s<=n) k_s exp(sum_(t<s) k_s) \ &= 1 + sum_(s<=n) k_s exp(sum_(t<s) k_s), $
as desired.

*Continuous case*. Given $g_0, k(t), p(t) >= 0$ and $ phi(t) <= g_0 + integral_0^t k(s) phi(s) + p(s) dif s, $ we need to prove $ phi(t) <= (g_0 + integral_0^t p(s) dif s) exp(integral_0^t k(s) dif s). $ Similarly, let $q(t) = g_0 + integral_0^t p(s) dif s$ be a non-decreasing function. We need to prove $ phi(t) <= q(t) + integral_0^t k(s) phi(s) dif s ==> phi(t) <= q(t) exp(integral_0^t k(s) dif s). $
For this let $ psi(t) = exp(-integral_0^t k(s) dif s)integral_0^t k(s) phi(s) dif s. $ Differentiating by the fundamental theorem of calculus, we get $ psi'(t) &= exp(-integral_0^t k(s) dif s) k(t) phi(t) - k(t) exp(-integral_0^t k(s) dif s)integral_0^t k(s) phi(s) dif s\
&= exp(-integral_0^t k(s) dif s) k(t) phi(t) - k(t) psi(t) \
&<= exp(-integral_0^t k(s) dif s) k(t) (q(t) + integral_0^t k(s) phi(s) dif s) - k(t) psi(t)\
&= exp(-integral_0^t k(s) dif s) k(t) q(t). $
Since $psi(0) = 0$, we can integrate this back.
$ psi(t) &<= integral_0^t exp(-integral_0^s k(r) dif r) k(s) q(s) dif s $
Substituting the definition of $psi$, and since $q$ is non-decreasing,
$ integral_0^t k(s) phi(s) dif s &<=  integral_0^t k(s) q(s) exp(integral_s^t k(r) dif r) dif s \
&<= q(t) integral_0^t k(s) exp(integral_s^t k(r) dif r) dif s  \
&= q(t) [exp(integral_0^t k(r) dif r) - 1] . $
Where the last equality is because the derivative of $F(s) = -exp(integral_s^t k(r) dif r)$ is exactly $F'(s) = k(s) exp(integral_s^t k(r) dif r)$, so the integral is evaluated by the fundamental theorem of calculus to be $F(t) - F(0)$.
Therefore we have
$ phi(t) &<= q(t) + integral_0^t k(s) phi(s) dif s\
&<= q(t) + q(t) [exp(integral_0^t k(r) dif r) - 1] = q(t) exp(integral_0^t k(r) dif r). $
This proves the result.

== Problem 2

For consistency, we suppose $t_0 = 0, t_1 = h$. It suffices to assume $u_0 = y_0 = 0$, and we analyze $tau(h) = |y_1 - u_1|/h$. Now we expand $f(t,y) = A + B t + C y + dots.c$. $ u_1 = h/6 [(4+h) A + 2f(h, u_1)] = 2/3 A h + 1/6 A h^2 + 1/3 h f(h, u_1) $
Solving this as a formal power series gives
$ u_1 &= A h + 1/6(A + 2 B - 2 A C) h^2 + O(epsilon^3) \
y_1 &= A h + 1/2(B + A C) h^2 + O(epsilon^3). $
Therefore the method is consistent, and the order of local truncation error is 1.

We test this method with two ODEs,
$ y' &= t + y \ y' &= -cos(y). $
We evaluate them at $t in [0,1]$. They both display a clean order 1 convergence:
#align(center, table(
    columns: (auto, auto, auto),
    [$h$], [$delta_1$], [$delta_2$],
    [$10^(-0)$], [$2.183 times 10^(-1)$], [$2.049 times 10^(-3)$],
    [$10^(-1)$], [$2.813 times 10^(-2)$], [$1.010 times 10^(-3)$],
    [$10^(-2)$], [$2.859 times 10^(-3)$], [$1.158 times 10^(-4)$],
    [$10^(-3)$], [$2.863 times 10^(-4)$], [$1.173 times 10^(-5)$],
    [$10^(-4)$], [$2.864 times 10^(-5)$], [$1.175 times 10^(-6)$],
    [$10^(-5)$], [$2.864 times 10^(-6)$], [$1.175 times 10^(-7)$],
    [$10^(-6)$], [$2.864 times 10^(-7)$], [$1.175 times 10^(-8)$],
))
Here $delta_1$ is the maximum error for the first ODE, and $delta_2$ that of the second.

== Problem 3

Without loss of generality, we assume $u_0 = y_0 = 0$. Expanding $f(t,y) = A + B t + C y + D t^2 + E t y + F y^2 + O(epsilon^3)$, we get
$ u_1 = A h  + (B + A C)/2 h^2 + (2D + 2 A E + B C + 2 A^2 F + C^2 A)/ 3 h^3 + O(h^4). $
We can also solve for $y = X t + Y t^2 + Z t^3 + O(t^4)$ to get $ y_1 =A h + (B + A C)/2 h^2 + (2 D + B C + A C^2 + 2 E A) / 6 h^3 + O(h^4).  $
Therefore the method is consistent with order of local truncation error 2. We just need to verify whether it is zero-stable. When $f, (diff f)/(diff t), (diff f)/(diff y)$ are all Lipschitz continuous in $y$, zero-stability is satisfied. Otherwise we can construct $f(t,y) = y^2 sin(1/y^2)$ with divergent (but defined) derivatives.

== Problem 4

Specializing to the Cauchy problem, we have
$ u_(n+2) - (1+alpha) u_(n+1) + alpha u_n = h/12 [(5 + alpha) lambda u_(n+2) + 8(1−alpha)lambda u_(n+1)−(1 + 5α)lambda u_n]. $
The characteristic equation is
$ (1 - (5+alpha)/12 h lambda) x^2 - (1+alpha + (2-2alpha)/3 h lambda) x + (alpha - (1+5alpha)/12 h lambda) = 0. $
We need $|x_1| , |x_2| < 1$. This is equivalent to $ |x_1|^2 + |x_2|^2 &< 2 \ |x_1|^2 + |x_2|^2 &< |x_1 x_2|^2 + 1. $ By Vieta's theorem we can simplify this to $ |b|^2 + |Delta| - 2|a|^2 < 2|c|^2 < 2 |a|^2. $
The latter simplifies to (when $|alpha|<1$) the region outside a circle with radius $1/2$ and center $(0,5/2)$. The former is a complicated shape:
#image("Hw9-Fig4.png")
Here all the regions contain at least an interval $(-epsilon, 0)$. For $alpha$ smaller than $1/7$, the region extends to positive and negative infinity. Otherwise the region is finite.

== Appendix: Code
Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

#raw(read("Hw9.py"), block: true, lang:"python")
