***************
Mathematics
***************

Numerical with Julia
====================

Contents:

- :ref:`machine-rep`
- :ref:`linear-alg`


.. _machine-rep:

Machine Representation of Numbers
--------------------------------------
``Julia`` supports multiple float-point types, including ``Float16`` (half), ``Float32`` (single) and  ``Float64`` (double). Half-precision floating-point numbers are only for storage format, when ``Float16`` type is involved in computation, it will be automatically promoted into ``Float32``. 

.. code-block:: julia
    :linenos:

    julia> sizeof(1.0)
    8
    julia> sizeof(float16(1.0))
    2
    julia> sizeof(float16(1.0)*2)
    4

``NaN`` and ``Inf`` (as well as ``-Inf``) are special float-points(IEC559), and can be cast into all float-point types also be used in arithmetic operations.

.. code-block:: julia
    :linenos:

    julia> Inf * Inf + Inf 
    Inf
    julia> Inf + (-Inf)
    NaN

Epsilon function ``eps`` gives machine accuracy. For single and double accuracy, epsilon will be ``float32(2.0^-23)`` and ``2.0^-53`` relatively. Notice ``eps`` also can be used on ``Float16``, the result will be also a half precision number.

.. code-block:: julia
    :linenos:

    julia> eps(Float16)
    float16(0.00097656)

    julia> eps(Float32)
    1.1920929f-7

    julia> eps(Float64)
    2.220446049250313e-16

The default rounding mode is ``RoundNearest``, to change the mode, we can use ``with_rounding``. In the following example, ``1.2000000000000001`` cannot be represented, thus ``Julia`` rounds it to the nearest representable float-point number.

.. code-block:: julia
    :linenos:

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
Linear algebra functions in ``Julia`` are largely implemented by calling functions from LAPACK. Sparse factorizations call functions from SuiteSparse.

- Matrices factorization
- Eigenvalues









Symbolic Computation
===========================