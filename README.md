# Numerical Methods in C

A collection of numerical methods implemented in **C** for **MATH 261** at Lake Superior State University.

The repository focuses on several substantial numerical-computing projects, including **PLU decomposition with partial pivoting, the Conjugate Gradient method, adaptive quadrature, and numerical evaluation of the normal cumulative distribution function**.

The projects emphasize implementing numerical algorithms directly, working with matrix operations and iterative methods, and comparing different approaches for accuracy and computational efficiency.

## Featured Projects

### PLU Decomposition

Implementation of matrix factorization using **partial pivoting** to decompose a matrix into permutation, lower-triangular, and upper-triangular components.

The project includes:

- PLU decomposition
- Partial pivoting
- Forward substitution
- Backward substitution
- Solution of linear systems
- Matrix operations
- Verification of computed solutions

The main numerical routines were implemented within a course-provided code framework.

```text
linear_algebra/plu_decomposition/
```

---

### Conjugate Gradient Method

Implementation of the **Conjugate Gradient method** for solving a large linear system arising from a heat/diffusion problem.

Rather than explicitly constructing and multiplying a large dense matrix, the implementation takes advantage of the known matrix structure to perform the required matrix-vector operations efficiently.

The project includes:

- Conjugate Gradient iteration
- Vector dot products
- Residual calculations
- Structured matrix-vector multiplication
- Execution-time measurement
- Comparison with direct linear-system solution methods

The implementation was tested on systems with sizes reaching approximately **100,000 unknowns**, demonstrating the scalability advantage of the iterative approach.

```text
linear_algebra/conjugate_gradient/
```

---

### Numerical Evaluation of the Normal CDF

A numerical approximation of the **normal cumulative distribution function** built using several numerical algorithms implemented in C.

The project combines:

- Newton's method for square-root approximation
- Numerical approximation of π using adaptive Simpson integration
- Numerical approximation of the exponential function
- Gaussian quadrature
- Adaptive numerical integration

These components are combined to numerically evaluate the normal distribution rather than relying entirely on standard-library implementations.

```text
numerical_integration/normal_cdf/
```

---

### Adaptive Gauss-Kronrod Quadrature

Implementation of an adaptive **7-15 Gauss-Kronrod quadrature method**.

The algorithm compares Gauss and Kronrod estimates to determine the local integration error and recursively subdivides intervals when additional accuracy is required.

The implementation is used to numerically approximate the natural logarithm through integration.

```text
numerical_integration/gauss_kronrod/
```

## Repository Structure

```text
Numerical_Methods_C/
│
├── README.md
│
├── linear_algebra/
│   │
│   ├── plu_decomposition/
│   │   ├── main.c
│   │   ├── dylanlib.h
│   │   ├── linalg.h
│   │   └── hardcode_matrix.h
│   │
│   └── conjugate_gradient/
│       ├── conjugate_gradient.c
│       ├── linalg.h
│       ├── heat.h
│       └── hires_timer.c
│
└── numerical_integration/
    │
    ├── normal_cdf/
    │   └── normal_cdf.c
    │
    └── gauss_kronrod/
        └── gauss_kronrod.c
```

## Numerical Methods Demonstrated

### Linear Algebra

- Matrix and vector operations
- PLU decomposition
- Partial pivoting
- Forward substitution
- Backward substitution
- Direct solution of linear systems
- Conjugate Gradient iteration
- Structured matrix-vector multiplication

### Numerical Integration

- Adaptive Simpson's rule
- Gaussian quadrature
- Gauss-Kronrod quadrature
- Adaptive error estimation

### Nonlinear Numerical Methods

- Newton's method
- Numerical approximation of elementary functions
- Error-controlled iterative algorithms

## Efficiency and Numerical Computing

Several projects in this repository explore not only whether a numerical method produces the correct result, but also how the choice of algorithm affects computational cost.

The Conjugate Gradient project is a particularly important example. For the structured linear system used in the heat problem, storing and multiplying the entire matrix would be unnecessarily expensive. Instead, the implementation computes the effect of the matrix directly from its known structure.

This allows substantially larger systems to be solved while using significantly less memory than a conventional dense-matrix representation.

Similarly, the adaptive integration projects evaluate local error estimates and refine only the regions that require additional computation.

## Technologies

- C
- Numerical Linear Algebra
- Iterative Methods
- Matrix Algorithms
- Numerical Integration
- Error Estimation
- Algorithm Analysis
- Dynamic Memory
- Performance Timing

## Course Context

These projects were completed for **MATH 261** at **Lake Superior State University**.

The repository contains a curated selection of the most substantial programming work from the course rather than every individual homework assignment.

Some projects were completed using instructor-provided frameworks, headers, or supporting code. Student-implemented numerical algorithms are retained alongside those supporting files where necessary to build or understand the programs.
