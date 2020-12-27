Erlang Native Matrix Functions
============================

This is an native Erlang linear algebra library (linalg)

# Installation
-----

Then, in the top directory, compile using rebar

> rebar compile

# Functions
-----

Creation and Description
 - shape(m)
 - row(m,i); col(m,j); cell(m,i,j)
 - zeros(i); zeros(i,j)
 - ones(i); ones(i,j)
 - sequential(i); sequential(i,j)
 - random(i); random(i,j)
 - eye(i); eye(i,j)
 - indentity(v)
 - diag(v); diag(m)

Matrix and Vector
 - dot(a,b)
 - inner(a,b)
 - outer(a,b)
 - transpose(m)
 - matmul(m,n)

Decomposition
 - qr(m) (for square matrix)
 - roots(v) (upto third order)
 # svd (under development)

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
 - log(m)

Matrix Reduction
 - sum, norm

# Usage
-----

These are a few examples. See the tests directory for more examples.

	Erlang R16B03 (erts-5.10.4) [source] [64-bit] 

	Eshell V5.10.4  (abort with ^G)

	1> linalg:dot([1.0,2.0],[3.0,4.0]).

	1> linalg:transpose([[1.0,2.0],[3.0,4.0]]).

	1> linalg:matmul([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).

    1> linalg:mul([1,2,3],[4,5,6]).
    
    1> linalg:add(10,[[1,2,3],[4,5,6]]).


License
-----
This project is licensed under the terms of the GNU General Public License v3.0	
