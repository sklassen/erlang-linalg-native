lapack - NIFs for Erlang
============================

This is an erlang matrix library 

Installation
-----

Then, in the top directory, compile using rebar

	> rebar compile

Usage
-----

Only two functions are supported as yet.

	Erlang R16B03 (erts-5.10.4) [source] [64-bit] 

	Eshell V5.10.4  (abort with ^G)
	1> matrix:t([[1.0,2.0],[3.0,4.0]]).
	1> matrix:mmultiply([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).



