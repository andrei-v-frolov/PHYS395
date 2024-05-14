# PHYS395 - Computational Physics

Simon Fraser University, Spring 2024

### Course description:

The course covers advanced numerical methods for scientific computing and provides introduction to programming in modern High Performance Computing (HPC) environment. Topics include:

- Representation of functions (cardinal vs. spectral basis, Fourier transform, orthogonal polynomials)
- Linear algebra (solving linear equations, least square fits, Cholesky and singular value decompositions)
- Statistics primer (generating random numbers, PDF and CDF estimators, MCMC and parameter likelihood)
- Root finding and optimization (bracketing and bisection, Newton's method, steepest descent and Levenberg-Marquardt algorithm)
- Neural networks (network topologies, network training, optimizing learning rate, PyTorch examples)
- Ordinary differential equations (integration methods, initial vs. boundary value problems)
- Hyperbolic partial differential equations (solution methods, numerical stability, wave equation)
- Parabolic PDEs (heat diffusion equation, numerical stability, spectral methods)
- Elliptic PDEs (boundary value problem revisited, Laplace equation, non-linear BVPs)
- Optimizing for performance; GPU acceleration and Fast Fourier Transforms revisited
- Scripting in Python (automating repeated tasks, making publication-quality plots)
- Going parallel on shared memory and MPI architectures (if time allows)

Programming environment is Python, with some exposure to compilers like Fortran later in the course. Homework and final are coding; you are expected to produce a working code that compiles, runs, and finds accurate numerical solution to the problem assigned. Bringing your own laptop is encouraged, but no technical support for Windows will be provided.

### More Info:

Lecture notes for the course are archived <a href="https://www.dropbox.com/scl/fo/jxe7wo5eedy1ygx6aihgz/AHitIFp8RjF55INcTFnO0ak?rlkey=9e3ozb1vuyfv7ivbydrzkhl30&dl=0">online</a>, and should serve as the code documentation of sorts, explaining what numerical techniques were used and why. For a quick view of what is covered in the course, take a look at the <a href="YY%20Gallery/README.md">gallery</a>. Achieving highest possible performance was not the main focus of the course, but most of the code was written with scalability in mind, and is fast enough to produce real-time animations and interactive demos on a laptop.
