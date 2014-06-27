***************
Mathematics
***************

Numerical with Julia
==========================

Contents:

- :ref:`machine-rep`
- :ref:`linear-alg`
- :ref:`optim`
- :ref:`approx`
- :ref:`quadrature`
- :ref:`ode`
- :ref:`pde`


.. _machine-rep:

Machine Representation of Numbers
--------------------------------------
``Julia`` supports multiple float-point types, including ``Float16`` (half), ``Float32`` (single) and  ``Float64`` (double). Half-precision floating-point numbers are only for storage format, when ``Float16`` type is involved in computation, it will be automatically promoted into ``Float32``. 

.. code-block:: julia


    julia> sizeof(1.0)
    8
    julia> sizeof(float16(1.0))
    2
    julia> sizeof(float16(1.0)*2)
    4

``NaN`` and ``Inf`` (as well as ``-Inf``) are special floating-points(IEC559), and can be cast into all floating-point types also be used in arithmetic operations.

.. code-block:: julia


    julia> Inf * Inf + Inf 
    Inf
    julia> Inf + (-Inf)
    NaN

Epsilon function ``eps`` gives machine accuracy. For single and double accuracy, epsilon will be ``float32(2.0^-23)`` and ``2.0^-53`` relatively. Notice ``eps`` also can be used on ``Float16``, the result will be also a half precision number.

.. code-block:: julia


    julia> eps(Float16)
    float16(0.00097656)

    julia> eps(Float32)
    1.1920929f-7

    julia> eps(Float64)
    2.220446049250313e-16

The default rounding mode is ``RoundNearest``, to change the mode, we can use ``with_rounding``. In the following example, ``1.2000000000000001`` cannot be represented, thus ``Julia`` rounds it to the nearest representable floating-point number.

.. code-block:: julia


    julia> 1.2000000000000001
    1.2000000000000002
    julia> with_rounding(Float64, RoundUp) do 
           1.1 + 0.1
           end
    1.2000000000000002
    julia> with_rounding(Float64, RoundDown) do 
           1.1 + 0.1
           end
    1.2

.. _linear-alg:

Linear Algebra
---------------
Linear algebra functions in ``Julia`` are largely implemented by calling functions from `LAPACK`_. Sparse factorizations call functions from SuiteSparse.

.. _LAPACK: http://www.netlib.org/lapack

Matrices factorization
^^^^^^^^^^^^^^^^^^^^^^

For factorization algorithms, ``Julia`` has already implemented common routines as ``lu``, ``qr`` and ``chol``. The routines call ``LAPACK`` subroutines for decomposition algorithm in general. For example, ``lufact!`` calls ``LAPACK.getrf!`` to update matrices in place.

.. code-block:: julia
    :emphasize-lines: 3

    function lufact!{T<:BlasFloat}(A::StridedMatrix{T};  pivot = true)
        !pivot && return generic_lufact!(A,pivot = pivot)
        lpt =  LAPACK.getrf!(A)
        return LU{T, typeof(A)}(lpt[1], lpt[2], lpt[3])

And ``LAPACK.getrf!`` calls ``getrf`` against the ``liblapack``.

.. code-block:: julia
    :emphasize-lines: 1-3

    ccall(($(string(getrf)), liblapack), Void, Ptr{BlasInt}, 
        Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
        Ptr{BlasInt}, &m, &n, A, &lda, ipiv, info)

Eigenvalues & Singular Value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The most common routines ``eig`` and ``svd`` call ``LAPACK.geevx!`` and ``LAPACK.gesdd!`` relatively. 

``svds`` does not exist but can be found with external package -- ``IterativeSolvers.jl``.

Solve ``Ax == b``
^^^^^^^^^^^^^^^^^
Solving linear system involves `direct method`_ and `iterative method`_.

``Julia`` has builtin function ``\(A,B)`` for solving linear system, which calls ``UMFPACK`` solver.

For iterative methods, there are plenty of routines(preconditioned) to choose from, such as ``CG``, ``GMRES``, ``SOR``, ``SSOR``, ``Lanczos``.

The above routines can be found through ``IterativeSolvers.jl`` for ``0.2.1-``. In this repository, we store our ``0.3.0+`` compatible routines in ``src`` directory, for example, following routines are stored in ``src/linalg/iterative``.  

.. literalinclude:: ../../src/linalg/iterative/cg.jl
   :language: julia
   :lines: 3-5

.. literalinclude:: ../../src/linalg/iterative/gmres.jl
   :language: julia
   :lines: 3-5

.. _direct method: http://en.wikipedia.org/wiki/direct_method

.. _iterative method: http://en.wikipedia.org/wiki/iterative_method


.. _optim:

Optimization
-------------
``JuliaOpt`` has built pure ``Julia`` package ``Optim`` for general optimization algorithms. Another package ``JuMP`` involves linear- and quadratic-constrained problems, which supports `Ipopt`_ through wrapper ``Ipopt.jl`` and many other solvers -- ``Cbc``, ``MOSEK``, ``Gurobi``, etc.

.. _Ipopt: https://projects.coin-or.org/Ipopt

.. _approx:

Approximation & Differentiation
--------------------------------




.. _quadrature:

Quadrature
--------------
``Julia`` has builtin support for one-dimensional integration as ``quadgk``, using adaptive *Gauss-Kronrod* algorithm. 

.. code-block:: julia
    :emphasize-lines: 1-4

    function quadgk{T<:FloatingPoint}(f, a::Complex{T},
                                      b::Complex{T}, c::Complex{T}...,
                                      abstol=zero(T), reltol=sqrt(eps(T))
                                      maxevals=10^7, order=7,norm=vecnorm)
        do_quadgk(f, [a,b,c...], order, T, abstol, reltol, maxevals, norm)
    end

For multi-dimensional case(including one dimensional), external package ``Cubature.jl`` can handle integrand over multi-dimensional box. 

In ``/src/integr``, we made more quadrature functions for use.


.. _ode:

Ordinary Differential Equation
-------------------------------

.. _pde:

Partial Differential Equation
------------------------------

Symbolic Computation
===========================