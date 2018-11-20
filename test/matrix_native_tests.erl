-module(matrix_native_tests). 
-import(matrix_native,[transpose/1,multiply/2]).
-include_lib("eunit/include/eunit.hrl").

transpose_test() ->
	transpose([[1.0,2.0],[3.0,4.0]])==[[1.0,3.0],[2.0,4.0]].

multiply_test()->
	multiply([[1.0,2.0],[3.0,4.0]],[[1.0,3.0],[2.0,4.0]])==[[5.0,11.0],[11.0,25.0]].

