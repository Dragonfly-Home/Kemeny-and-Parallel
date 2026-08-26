# Kemeny-and-Parallel

MATLAB code accompanying the paper "On a Decomposition Formula for
Kemeny's Constant via Stochastic Complement."

## Third-Party Code

The implementations of the direct and recursive baseline methods in

- `src/kemenydirect.m`
- `src/recursivekemeny.m`
- `src/recursivekemenydirect.m`

are reproduced from the
[Kemeny-and-Conquer](https://github.com/Cirdans-Home/Kemeny-and-Conquer)
repository by Dario Andrea Bini, Fabio Durastante, Sooyeong Kim, and
Beatrice Meini.

These files remain subject to the BSD 3-Clause License. The original
copyright notice and license conditions are retained in
`LICENSE-Kemeny-and-Conquer`. See `THIRD_PARTY_NOTICES.md` for details.

## Data

The test matrices are obtained from the SuiteSparse Matrix Collection.
Several matrix files are also distributed with the Kemeny-and-Conquer
repository.

## License

Original code developed for Kemeny-and-Parallel is licensed under the
MIT License in `LICENSE`. Third-party code identified above is licensed
under the BSD 3-Clause License in `LICENSE-Kemeny-and-Conquer`.
