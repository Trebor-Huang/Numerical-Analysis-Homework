#set text(font: "CMU Serif")
= Homework 7

Time to try out some Typst!

== Problem 1

#set enum(numbering: "(a)")
+ We need to calculate $ F'(x) = 1 - ((f'(x))^2 - f(x)f''(x))/(f'(x)^2) = (f(x)f''(x))/(f'(x)^2). $ By Taylor expansion, we have $ f(x) &= (x-alpha)^m (r + o(1)) \
f'(x) &= (x-alpha)^(m-1) (m r + o(1)) \
f''(x) &= (x-alpha)^(m-2) (m(m-1)r + o(1)). $ Note that these are independent Taylor expansions, i.e. the latter does not come from a derivative of the former, since $o(1)$ is generally not differentiable. Therefore we can substitute and expand to get $ F'(x) &= ([(x-alpha)^m (r + o(1))][(x-alpha)^(m-2) (m(m-1)r + o(1))]) / ([(x-alpha)^(m-1) (m r + o(1))]^2) \
&= ([r + o(1)] [m(m-1)r + o(1)]) / ([m r + o(1)]^2). $ Taking the limit $x -> alpha$, we get $F'(x) -> (r m (m-1) r) / ((m r)^2) = (m-1)/m$, as desired.

+ We only need to confirm that the modified $F$ satisfies $F'(alpha) = 0$. (This ensures that quadratic convergence is _achieved_, but in special cases it may be _surpassed_ too.) Now $F(x) = x - m(x - G(x))$, where $G$ is the unmodified version. Therefore by linearity of derivatives, $ F'(alpha) = 1 - m(1 - G'(alpha)) = 1 - m(1 - (m-1)/m) = 0. $

== Problem 2

#set enum(numbering: "(i)")
+ We are iterating $ F(x) = x - (f(x)^2) / (f(x+f(x)) - f(x)). $ Let $alpha$ be the root in question. We need to prove $F(alpha) = alpha$ and $F'(alpha) = 0$. For the first claim, we use the L'Hôpital rule to get $ lim_(x -> alpha) F(x) &= alpha - lim_(x -> alpha) (2f(x)f'(x)) / (f'(x+f(x))(1 + f'(x)) - f'(x))\ &= alpha - (0) / (f'(alpha)(1 + f'(alpha)) - f'(alpha)) = alpha $ For the second equation, differentiating we have $ F'(x) &= 1 - (2f(x)f'(x)) / (f(x+f(x)) - f(x)) + (f(x)^2[f'(x+f(x))(1+f'(x)) - f'(x)] ) / ([f(x + f(x)) - f(x)]^2). $ Using the L'Hôpital rule, we have $ lim_(x -> alpha) (2f(x)f'(x)) / (f(x+f(x)) - f(x)) &= lim_(x -> alpha) (2(f(x)f''(x) + f'(x)f'(x))) / (f'(x+f(x))(1 + f'(x)) - f'(x))\ &= (2f'(alpha)^2) / (f'(alpha)(1 + f'(alpha)) - f'(alpha)) \ &= 2. $ For the other term, we have $ &lim_(x -> alpha) (f(x)^2[f'(x+f(x))(1+f'(x)) - f'(x)] ) / ([f(x + f(x)) - f(x)]^2) \ =& [f'(alpha)(1+f'(alpha)) - f'(alpha)] [lim_(x -> alpha) (f(x)) / (f(x + f(x)) - f(x))]^2 \ =& f'(alpha)^2 [lim_(x -> alpha) (f'(x)) / (f'(x + f(x))(1 + f'(x)) - f'(x))]^2 \ =& f'(alpha)^2 [ (f'(alpha)) / (f'(alpha)^2) ]^2 = 1. $ Putting everything together, $lim_(x -> alpha) F'(x) = 1 - 2 + 1 = 0$. This proves the claim.

+ There is cancellation at $f(x+f(x)) - f(x)$. So we have to stop when this number becomes zero (i.e. smaller than machine precision). It gives the following results: #figure(table(columns: (auto, auto),
[$x_0=0$], [$x_0=3$],
[$9.00814 times 10^(-2)$], [$1.83875 times 10^(-3)$],
[$1.12561 times 10^(-3)$], [$3.06533 times 10^(-7)$],
[$1.95964 times 10^(-7)$], [$8.88178 times 10^(-15)$],
[$5.99520 times 10^(-15)$], [$4.44089 times 10^(-16)$],
[0], [0],), caption: [Errors for each iteration]) The result does suggest quadratic convergence.

+ We can consider $ F_1(x) = arcsin(e^(-x)) $ (which restricts $x in [-1, 1]$, which is good because there is only one root in this region), and $ F_2(x) = e^(-x) - sin(x) + x. $ For $F_1$, the derivative of $F'_1(alpha) approx -0.667$, suggesting linear convergence. $F_2$ is particularly interesting. The derivative is $F'_2(x) = -e^(-x) - cos(x) + 1$. For large $x$, the first factor is neglegible. So $F'_2(x_n)$ is approximately zero for odd roots, and approximately $2$ for even roots. Therefore it will never converge to even roots, while it converges linearly but very quickly for odd roots. #figure(image("Hw7-Fig1.png", width: 60%) + image("Hw7-Fig2.png", width: 60%)) Indeed they show linear convergence at drastically different speeds.

== Problem 3

We know the maximum is attained at $x=0$ but let's pretend we don't. Newton's method requires $f'(x)$ and $f''(x)$, which can be readily computed. Choosing an initial value of $x_0 = 1$, we see that it has linear convergence. The plot is just visually identical to a straight line so I will not present it here.

For an improvement, we can introduce a factor $m=3$, since $f'$ has a third-order root at the origin. This makes the convergence almost instant.

#pagebreak()

== Appendix: Code

Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

```python
import numpy as np
# Problem 2.

def f(x):
    return np.exp(-x) - np.sin(x)

def steffensen(f, x):
    while (u := f(x+f(x)) - f(x)) != 0:
        x = x - f(x) ** 2 / u
        yield x

print("\nSteffensen, Initial condition x=0:")
solution = steffensen(f, 0)
exact = 0.58853274398186107743245205
for s in solution:
    print("%.5e" % abs(s - exact))

print("\nSteffensen, Initial condition x=3:")
solution = steffensen(f, 3)
exact = 3.09636393241064611562584085
for s in solution:
    print("%.5e" % abs(s - exact))

print("\nIteration for F1:")
exact = 0.58853274398186107743245205
x = 0.5
for _ in range(20):
    x = np.arcsin(np.exp(-x))
    print("%.5e" % abs(x - exact))

print("\nIteration for F2:")
exact = 6.2850492733825866764618695
x = 6
for _ in range(10):
    x = np.exp(-x) - np.sin(x) + x
    print("%.5e" % abs(x - exact))

# Problem 3

def f(x):
    return np.exp(-x**2/2) - np.cos(2*x)/4
def df(x):
    return -x*np.exp(-x**2/2) + np.sin(2*x)/2
def ddf(x):
    return (x**2 - 1) * np.exp(-x**2/2) + np.cos(2*x)

print("\nNewton's method for maxima:")
x = 1
for _ in range(20):
    x = x - df(x) / ddf(x)
    print("%.5e" % x)

print("\nCorrected Newton's method for maxima:")
x = 1
for _ in range(5):
    x = x - 3 * df(x) / ddf(x)
    print("%.5e" % x)
```
