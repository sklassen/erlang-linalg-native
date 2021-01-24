-module(linalg_matrix_tests).
-import(linalg,[zeros/1,zeros/2,ones/1,ones/2,fill/2,fill/3]).
-import(linalg,[det/1,inv/1,transpose/1,dot/2,matmul/2]).
-include_lib("eunit/include/eunit.hrl").

dot_3_test() ->
	?assertEqual(32.0,dot([1.0,2.0,3.0],[4.0,5.0,6.0])).

transpose_test() ->
    [
	?assertEqual([[8.0]],transpose([[8.0]])),
	?assertEqual([[1.0,3.0],[2.0,4.0]],transpose([[1.0,2.0],[3.0,4.0]])),
	?assertEqual([[1.0,4.0],[2.0,5.0],[3.0,6.0]],transpose([[1.0,2.0,3.0],[4.0,5.0,6.0]]))
    ].

matmul_test()->
    [
    ?assertEqual([[6.0]],matmul([[2.0]], [[3.0]])),
    ?assertEqual([[5.0,11.0],[11.0,25.0]],matmul([[1.0,2.0],[3.0,4.0]], [[1.0,3.0],[2.0,4.0]])),
    ?assertEqual([[30, 36, 42]],matmul([[1,2,3]], [[1,2,3],[4,5,6],[7,8,9]])),
    ?assertEqual([[14], [32], [50]],matmul([[1,2,3],[4,5,6],[7,8,9]], [[1], [2], [3]])),
    ?_assertException(error, function_clause, matmul([[2.0], [3]], [[3.0]]))
    ].

det_test()->
    [
	?assertEqual(6,det([[6]])),
	?assertEqual(-2,det([[1,2],[3,4]])),
	?assertEqual(-306,det([[6,1,1], [4,-2,5], [2,8,7]]))
    ].


inv_test()->
    [
	?assertEqual([[0.125]],inv([[8]])),
	?assertEqual([[1.0,0.0],[0.0,0.5]],inv([[1,0],[0,2]])),
	?assertEqual([[-1.0,-1.0,2.0],[-1.0,0.0,1.0],[2.0,1.0,-2.0]],inv([[1,0,1],[0,2,1],[1,1,1]]))
    ].

zeros_test()->
    [
    ?assertEqual([],zeros(0)),
    ?assertEqual([0,0],zeros(2)),
    ?assertEqual([[]],zeros(0,0)),
    ?assertEqual([[]],zeros(1,0)),
    ?assertEqual([[]],zeros(0,2)),
    ?assertEqual([[0]],zeros(1,1)),
    ?assertEqual([[0,0,0],[0,0,0]],zeros(2,3))
    ].

ones_test()->
    [
    ?assertEqual([],ones(0)),
    ?assertEqual([1,1],ones(2)),
    ?assertEqual([[]],ones(0,0)),
    ?assertEqual([[]],ones(1,0)),
    ?assertEqual([[]],ones(0,2)),
    ?assertEqual([[1]],ones(1,1)),
    ?assertEqual([[1,1,1],[1,1,1]],ones(2,3))
    ].

fill_test()->
    [
    ?assertEqual(zeros(2),fill(2,0)),
    ?assertEqual(ones(3),fill(3,1)),
    ?assertEqual([2,2,2,2],fill(4,2)),
    ?assertEqual([[]],fill(0,0,8)),
    ?assertEqual([[]],fill(1,0,7)),
    ?assertEqual([[]],fill(0,2,1)),
    ?assertEqual(zeros(2,3),fill(2,3,0)),
    ?assertEqual(ones(3,4),fill(3,4,1)),
    ?assertEqual([[2],[2],[2],[2]],fill(4,1,2))
    ].
