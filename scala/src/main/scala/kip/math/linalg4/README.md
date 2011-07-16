
** Features **

- High performance complex number support; backed by primitive buffer, not heap allocated

- Blas and Lapack operations call native libraries using JNA

- Minimal dependencies

- Extensible typeclass based design

- Multiple matrix types (Dense, DenseRow/DenseColumn, Sparse, Hermitian, Triangular, Tridiagonal, Unitary)

- Choice of in-place operations, or natural math syntax


** Non-goals **

- Special support for small matrices (dimensions 2 or 3)

- No features beyond linear algebra (no plotting, general purpose implicits)
