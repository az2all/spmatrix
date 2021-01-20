<h1>Sparse Matrix Library in C++ Using Old Yale Format (<b>SMLC++</b>)</h1>.

<i><b>SMLC++</b></i> is a fully cross-platform library, which runs on Windows and Linux.

A sparse matrix is a matrix with large number of zero elements. An efficient way to store it
saving only non-zero elements is provided by the Old Yale format [1]. It reduces redundancy in
arrays by storing row and column information in a compact way.

Consider the <i>NxN</i> system of linear equations given by:<br>
<b>A</b> &middot; <b>x</b> = <b>b</b><br>
where the coefficient matrix <b>A</b> is sparse, large, symmetric, and positive definite.

<i><b>SMLC++</b></i> allows us to solve the system of equations using one of the following iterative methods:
- Jacobi,
- Gauss-Seidel,
- SOR (Successive Overrelaxation),
- Conjugate gradient.

main.cpp provides an example of library use.

<br>
<h3>References</h3><br>
[1] Eisenstat, S. C.; Gursky, M. C.; Schultz, M. H.; Sherman, A. H. <i>Yale
sparse matrix package I: the symmetric codes</i>, Int. J. Numer. Methods in Engin.,
August 1982, vol. 18, issue 8, pp. 1145-1151.
