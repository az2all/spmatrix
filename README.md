<h1>Sparse Matrix Library in C++ Using Old Yale Format (<b>SMLC++</b>)</h1>.

<b>SMLC++</b> is a fully cross-platform library.

A sparse matrix is a matrix with large number of zero elements. An efficient way to store it
saving only non-zero elements is provided by the Old Yale format. It reduces redundancy in
arrays by storing row and column information in a compact way.

Assume a system of linear equations given by:<br>
<b>A</b> &middot; <b>x</b> = <b>b</b><br>
where <b>A</b> is a large SPD sparse matrix.

<b>SMLC++</b> allows us to solve the system of equations using one of the following iterative methods:
- Jacobi,
- Gauss-Seidel,
- SOR (Successive Overrelaxation),
- Conjugate gradient.

main.cpp provides several examples of library use.

<h3>References</h3><br>
S. C. Eisenstat and M. C. Gursky and M. H. Schultz and A. H. Sherman, Yale
sparse matrix package I: the symmetric codes, Int. J. Numer. Methods in Engin.,
18 (1982), pp. 1145-1151.
