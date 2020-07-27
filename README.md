Erlang Native Matrix Functions
============================

This is an Erlang linear algebra library (linalg)

Installation
-----

Then, in the top directory, compile using rebar

	> rebar compile

Usage
-----

Only a few functions are supported at the moment, see the tests directory for examples.

	Erlang R16B03 (erts-5.10.4) [source] [64-bit] 

	Eshell V5.10.4  (abort with ^G)

	1> linalg_native:dot([1.0,2.0],[3.0,4.0]).

	1> linalg_native:transpose([[1.0,2.0],[3.0,4.0]]).

	1> linalg_native:matmul([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).



