<h1>Sparse Matrix Library in C++ Using Old Yale Format (<b>SMLC++</b>)</h1>.

<i><b>SMLC++</b></i> is a fully cross-platform library, which runs on Windows and Linux.

A sparse matrix is a matrix with large number of zero elements. An efficient way to store it
saving only non-zero elements is provided by the Old Yale format [1]. It reduces redundancy in
arrays by storing row and column information in a compact way.

Consider the <i>NxN</i> system of linear equations given by:<br>
<b>A x</b> = <b>b</b><br>
where the coefficient matrix <b>A</b> is sparse, large, symmetric, and positive definite.

<i><b>SMLC++</b></i> allows users to solve the system of linear equations using one of the following iterative methods:
- Jacobi,
- Gauss-Seidel,
- SOR (Successive Overrelaxation),
- Conjugate gradient.

main.cpp provides an example of library use.
<br><br>
<b>Example 1: Compute importance of graph nodes (PageRank) using power method</b>
![graph](https://user-images.githubusercontent.com/77605006/106374154-642cb480-6335-11eb-81e0-53267752c9d3.png)
<br>A reasonable solution is achieved after 5 iterations: v=[0.07 0.33 0.21 0.11 0.21 0.08].
<br><br>
<b>Example 2: System of linear equations Ax=b</b><br>
![graph2](https://user-images.githubusercontent.com/77605006/106374707-6d6c5000-633a-11eb-8bbc-3467edd36125.png)
<br>
CG method provides solution: <b>x</b>=[0.90 1.23 0.68 0.08].


<br>
<h3>References</h3><br>
[1] Eisenstat, S. C.; Gursky, M. C.; Schultz, M. H.; Sherman, A. H. <i>Yale
sparse matrix package I: the symmetric codes</i>, Int. J. Numer. Methods in Engin.,
August 1982, vol. 18, issue 8, pp. 1145-1151.<br>
[2] http://www.codedevelopment.net/algorithms.html
