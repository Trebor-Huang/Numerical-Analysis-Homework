#set text(font: "CMU Serif")
= Homework 11

== Problem 1

(i) The inverse of the matrix is $ A^(-1) = mat(375,-187;-376,375/2). $

We know that the matrix norms: $||A||_2$ is the square root of the largest eigenvalue of $A^*A$. $||A||_oo$ is the maximum absolute row sum. So we get the condition number $ kappa_oo(A) = ||A||_oo ||A^(-1)||_oo = 1502 times 563.5 = 846377. $ $ kappa_2(A) = ||A||_2 ||A^(-1)||_2 = sqrt((1983886335017 + 12676545 sqrt(24492423889))/8) approx 700000. $
This exercise is so cruel.

(ii) We have $||A#[]^(-1)||#[]^(-1) |v| <= | A v | <= ||A|| |v|$. Therefore we want $x$ to lean on the right side, and $delta x$ on the left. Thus, take $ delta b = vec(-1,1), x = vec(1,1). $ This gives $delta x = vec(-562,563.5), b = vec(749,1502)$, so $ ||delta b||/||b|| = 1/1502, quad ||delta x||/||x|| = 1127/2. $ They differ by a factor of $kappa_oo(A) = 846377$, as expected.

(iii) Just do the reverse: $ b = vec(-1,1), delta x = vec(1,1). $ This gives $x = vec(-562,563.5), delta b = vec(749,1502)$, so $ ||delta b||/||b|| = 1502, quad ||delta x||/||x|| = 2/1127. $

== Problem 2

(i) Since the first column has only two non-zero elements, either pivoting doesn't happen or it happens between the first two rows.

If pivoting doesn't happen, then $|a_(1,1)| >= |a_(2,1)|$. The elimination only affects one element, $a'_(2,2) = a_(2,2) - a_(1,2) a_(2,1)/a_(1,1)$. And by the previous inequality,
$ |a'_(2,2)| <= |a_(2,2)| + |a_(1,2)| |a_(2,1)|/|a_(1,1)| <= |a_(2,2)| + |a_(1,2)| <= 2 max|a_(i,j)|. $
Now by induction, the submatrix has a growth factor not exceeding 2. So $ max_(i,j>1) |u_(i,j)|<= 2 max_(i,j>1) |a'_(i,j)| = 2max { max_("except" a_(2,2))|a_(i,j)|, |a'_(2,2)| }. $

If pivoting does happen, there are two elements affected, and the submatrix remains a tridiagonal matrix. Thus induction still applies, and we can get a similar estimate.

(ii) We argue by induction that no pivoting is performed in this case. For the first column, since the matrix is column diagonally dominant, the pivot point is chosen to be $a_(1,1)$ itself. now we multiply the first row by $- a_(k,1)/a_(1,1)$, and add it to the $k$-th row. We still needs to prove that after removing the first row and column from the resulting matrix we still have a column diagonally dominant matrix. We originally have $ |a_(k,k)| >= sum_(j!=k) |a_(j,k)|. $ After the modification we have $a'_(i,j) = a_(i,j) - a_(i,1)/a_(1,1) a_(1,j)$ except for the first row which is not changed. So $ |a'_(k,k)| &= |a_(k,k) - a_(k,1)/a_(1,1) a_(1,k)| \ &>= |a_(k,k)| - |a_(k,1)|/|a_(1,1)| |a_(1,k)| \ &>= [sum_(j != k) |a_(j,k)| ] - |a_(k,1)|/|a_(1,1)| |a_(1,k)| $
And (note that we removed the first row) $ sum_(j != 1,k) |a'_(j,k)| &= sum_(j != 1,k) |a_(j,k) - a_(j,1)/a_(1,1) a_(1,k)| \ &<= sum_(j != 1,k) [ |a_(j,k)| + |a_(j,1)|/|a_(1,1)| |a_(1,k)| ] \ &= [ sum_(j != 1,k)  |a_(j,k)| ] + sum_(j != 1) [ |a_(j,1)|/|a_(1,1)| |a_(1,k)| ] - |a_(k,1)|/|a_(1,1)| |a_(1,k)| \ &<= [ sum_(j != 1,k)  |a_(j,k)| ] + (|a_(1,1)| - |a_(k,1)|)/|a_(1,1)| |a_(1,k)| \ &= [ sum_(j != k)  |a_(j,k)| ] - (|a_(k,1)|)/|a_(1,1)| |a_(1,k)|. $ Therefore the inequality is established.

To estimate the growth factor, we just need prove that $max_(i,j) |u_(i,j)|$ is attained on the diagonal, because it follows from the Gershgorin circle theorem that each one of the eigenvalues (i.e. the diagonal entries of $u$) is within a disk of center $a_(i,i)$ and radius $sum_(j!=i)|a_(j,i)|$, which is no more than $|a_(i,i)|$, and therefore $|u_(k,k)| <= 2 max_i|a_(i,i)|$. Since $A = L U$, and the entries of $L$ are no more than 1, if $max_(i,j) |u_(i,j)|$ is not attained on the diagonal then $L U$ cannot be diagonally dominant. This concludes the proof.


== Problem 3

(i) We have
#align(center, table(columns: (auto, auto),
[$n$], [$kappa(H_n)$],
[5], [$9.437 times 10^(5)$],
[6], [$2.907 times 10^(7)$],
[7], [$9.852 times 10^(8)$],
[8], [$3.387 times 10^(10)$],
[9], [$1.100 times 10^(12)$],
[10], [$3.536 times 10^(13)$],
[11], [$1.234 times 10^(15)$],
[12], [$4.016 times 10^(16)$],
[13], [$1.658 times 10^(18)$],
[14], [$5.488 times 10^(17)$],
[15], [$4.661 times 10^(17)$],
[16], [$1.239 times 10^(18)$],
[17], [$5.498 times 10^(17)$],
[18], [$3.673 times 10^(19)$],
[19], [$2.737 times 10^(19)$],
[20], [$9.131 times 10^(18)$],))
It grows exponentially until reaching around $10^19$. At this point the growth is curbed by machine precision.

(ii) Since the growth factor goes as $2^n$, we expect the error to grow at that rate too. We take the normal distribution around the origin. Averaging the error over 50 vectors, we get the result:
#align(center, image("Hw12-Fig1.png"))
This confirms our predictions.

== Appendix: Code
Code can also be found at #link("https://github.com/Trebor-Huang/Numerical-Analysis-Homework", underline[the GitHub repository]).

#raw(read("Hw12.py"), block: true, lang:"python")
