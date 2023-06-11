#set text(font: "CMU Serif")
= Homework 13

== Problem 1

Taking the singular value decomposition $A = U Sigma V^*$, we have $ A^+ &= (A^* A)^(-1) A^* = (V Sigma^* U^* U Sigma V^*)^(-1) (V Sigma^* U^*)\ &= V (Sigma^* Sigma)^(-1) Sigma^* U^* = V Sigma^+ U^*. $
Since $Sigma$ is a diagonal rectagular matrix, we can easily compute the pseudoinverse: Just transpose and invert elementwise. Therefore $||A^+||_2$ is its largest singular value, i.e. the reciprocal of $A$'s smallest singular value. On the other hand, $||A_1^(-1)||_2 = ||A_1^+||_2$ is the reciprocal of $A_1$'s smallest singular value. Note that what $A v$ does to the vector $v$ is it concatenates $A_1 v$ and $A_2 v$. So the norm can only increase, and therefore $A$ will always stretch the vectors more, i.e. have a larger minimal singular value. This proves the claim.

== Problem 2

We use an equidistant sampling, so the least-square method doesn't need any scaling to be consistent with the $L^2$ norm.

The optimal fit is
$ y = c_1 e^x + c_2 sin x + c_3 Gamma(x) $
where $c_1 = -0.10774242, c_2 = 0.00959252, c_3 = 1.28659552$. This produces the following plot:
#align(center, image("Hw13-Fig1.png", width:60%))

For $[0,1]$ the question doesn't even make sense because $1/x$ is simply not in $L^2[0,1]$. This can be fixed by demanding $1/x - c_3 Gamma(x)$ be in $L^2$, this determines $c_3 = 1$. The same procedure gives $c_1 = 0.6238602, c_2 = -1.81576689$.

== Problem 3

Implementation details commented in code.

Notable things: Usually $"sgn"(0) = 0$, but in this particular application, we want $"sgn"(0) = 1$ (or $-1$, doesn't matter) because we are choosing which side we want to reflect to. Even if the vector is right in the middle, we would need to choose one side.

The six methods give relative errors of $x$:
- Cholesky: $5.09 times 10^(-2)$
- QR factorization methods:
    - Modified Gram--Schmidt: $1.25 times 10^(-3)$
    - Householder transform: $4.33 times 10^(-10)$
    - Householder with pivoting (built-in): $1.01 times 10^(-9)$
- Pseudoinverse: $3.73 times 10^(-10)$
- Singular value decomposition: $3.05 times 10^(-10)$.
It is numerically checked that there are no computation errors, for example $Q R = A$ indeed holds for the MGS algorithm. So the differences are purely due to numeric error.

It is seen that directly solving the normal equations behaves very badly numerically. Modified Gram--Schmidt also loses some precision on the small components.

== Problem 4

(I'm sorry, but I'm a bit too tight on final deadlines. I promise I'll do it afterwards.)

== Appendix: Code
Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

#raw(read("Hw13.py"), block: true, lang:"python")
