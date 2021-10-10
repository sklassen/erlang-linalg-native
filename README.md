Erlang Native Matrix Functions
============================

This is an native Erlang linear algebra library (linalg). 
The function signatures and results strive to match numpy and scipy linalg.

The are also erlang-linalg NIF libraries avalible that match in function 
names, signature and results. This is the native version is slower than
those NIFs. This code was written to be readable and (maybe) teachable rather
than fast.

Erlang's lists:nth indexes start at 1; and linalg's Row, Col, Cell, ArgMax and ArgMin 
function use this convention. This differs from python et al at start at 0.
Like the numpy and scipy, matrix input is alway expected to be Row Major. 

Suggestions, enhancements and pull-requests welcomed. 

# Installation

There are no dependancies. In the top directory, compile using rebar.

> rebar compile

To include it as a rebar.config dependancy, add the line below.
```
{deps,[
    {linalg, ".*", {git, "https://github.com/sklassen/erlang-linalg-native.git"}}
]}
```

# Functions

Creation and Description
 - shape(m)
 - row(i,m); col(j,m); cell(i,j,m); set_cell(i,j,value,m)
 - zeros(i); zeros(i,j)
 - ones(i); ones(i,j)
 - sequential(i); sequential(i,j)
 - random(i); random(i,j)
 - fill(i,value); fill(i,j,value)
 - eye(i); eye(i,j)
 - identity(v)
 - diag(v); diag(m)

Matrix and Vector
 - dot(a,b)
 - inner(a,b)
 - outer(a,b)
 - transpose(m)
 - matmul(m,n)

Decomposition
 - lu(m)
 - qr(m) 
 - cholesky(m)
 - roots(v) (upto third order)

Norms and others
 - norm(a)
 - det(m)

Solves and inversion
 - inv(m)
 - solve(m)

Basic Scalar, Vector and Matrix Arithmetric 
 - add(n,m)
 - sub(n,m)
 - mul(n,m)
 - divide(n,m)
 - pow(n,m)
 - sqrt(m)
 - abs(m)
 - log(m)

Reductions
 - sum(m)
 - sumsq(m)
 - norm(m)
 - mean(m)
 - median(m)
 - std(m)
 - var(m)
 - cov(m)
 - min(m)
 - max(m)
 - argmin(m)
 - argmax(m)

# Usage

Below are a few examples. See the tests directory for more examples.
```
Erlang R16B03 (erts-5.10.4) [source] [64-bit] 

Eshell V5.10.4  (abort with ^G)

1> linalg:dot([1.0,2.0],[3.0,4.0]).

1> linalg:transpose([[1.0,2.0],[3.0,4.0]]).

1> linalg:matmul([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).

1> linalg:mul([1,2,3],[4,5,6]).
    
1> linalg:add(10,[[1,2,3],[4,5,6]]).
```

# TODO
 - adjoint
 - eig
 - eigh
 - eigvals
 - eigvalsh
 - lstsq 
 - svd
 - pinv

# License

This project is licensed under the terms of the GNU General Public License v3.0	
