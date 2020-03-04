-module(linalg_native_tests). 
-import(linalg_native,[det/1,inv/1,transpose/1,dot/2,matmul/2]).
-include_lib("eunit/include/eunit.hrl").

dot_3_test() ->
	?assert(dot([1.0,2.0,3.0],[4.0,5.0,6.0])==32).

transpose_1_test() ->
	?assert(transpose([[8.0]])=:=[[8.0]]).

transpose_2_test() ->
	?assert(transpose([[1.0,2.0],[3.0,4.0]])==[[1.0,3.0],[2.0,4.0]]).

transpose_3_test() ->
	?assert(transpose([[1.0,2.0,3.0],[4.0,5.0,6.0]])==[[1.0,4.0],[2.0,5.0],[3.0,6.0]]).

matmul_1x1_test()->
	?assert(matmul([[2.0]],[[3.0]])==[[6.0]]).

matmul_2x2_test()->
	?assert(matmul([[1.0,2.0],[3.0,4.0]],[[1.0,3.0],[2.0,4.0]])==[[5.0,11.0],[11.0,25.0]]).


det_1_test()->
	?assert(det([[6]])==6).

det_2_test()->
	?assert(det([[1,2],[3,4]])==-2.0).

det_3_test()->
	?assert(det([[6,1,1], [4,-2,5], [2,8,7]])==-306.0).


inv_1_test()->
	?assert(inv([[8]])==[[0.125]]).

inv_2_test()->
	?assert(inv([[1,0],[0,2]])==[[1.0,0.0],[0.0,0.5]]).

inv_3_test()->
	?assert(inv([[1,0,1],[0,2,1],[1,1,1]])==[[-1.0,-1.0,2.0],[-1.0,0.0,1.0],[2.0,1.0,-2.0]]).

