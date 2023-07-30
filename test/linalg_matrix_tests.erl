-module(linalg_matrix_tests).
-import(linalg, [cell/3, set_cell/4, set_row/3, set_col/3]).
-import(linalg, [zeros/1, zeros/2, ones/1, ones/2, fill/2, fill/3]).
-import(linalg, [inv/1, transpose/1, dot/2, matmul/2]).
-include_lib("eunit/include/eunit.hrl").

cell_test() ->
    Matrix = [[1, 2, 3], [4, 5, 6]],
    [
        ?assertEqual(cell(1, 1, Matrix), 1),
        ?assertEqual(cell(2, 3, Matrix), 6),
        ?assertException(error, function_clause, cell(10, 2, Matrix))
    ].

set_cell_test() ->
    Matrix = [[1, 2, 3], [4, 5, 6]],
    [
        ?assertEqual(set_cell(1, 1, 2, Matrix), [[2, 2, 3], [4, 5, 6]]),
        ?assertEqual(set_cell(2, 2, 4, Matrix), [[1, 2, 3], [4, 4, 6]]),
        ?assertException(error, function_clause, set_cell(10, 2, 3, Matrix))
    ].

set_row_test() ->
    Matrix = [[1, 2, 3], [4, 5, 6]],
    [
     ?assertEqual(set_row(1, [7, 8, 9], Matrix), [[7, 8, 9], [4, 5, 6]])
    ].

set_col_test() ->
    Matrix = [[1, 2, 3], [4, 5, 6]],
    [
     ?assertEqual(set_col(1, [7, 8], Matrix), [[7, 2, 3], [8, 5, 6]])
    ].

dot_0_test() ->
    ?assertEqual(12, dot(3,4)).
dot_1_test() ->
    ?assertEqual(32.0, dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])).
dot_2_test() ->
    ?assertEqual([[4, 1],[2, 2]], dot([[1, 0], [0, 1]],[[4, 1], [2, 2]])).

dot_3_test() ->
    ?assertEqual([[8, 2], [4, 4]], dot([[4, 1], [2, 2]],2)).
dot_4_test() ->
    ?assertEqual([[8, 2], [4, 4]], dot(2,[[4, 1], [2, 2]])).

dot_5_test() ->
    ?assertEqual([11, 10], dot([[4, 1], [2, 2]],[2,3])).
dot_6_test() ->
    ?assertEqual([14, 8], dot([2,3],[[4, 1], [2, 2]])).

dot_7_test() ->
    ?assertEqual([[0.22000000000000003,0.2475,0.0825]], dot([[0.5, 0.25, 0.25]],[[0.4, 0.03, 0.05], [0.03, 0.9, 0.03], [0.05, 0.03, 0.2]])).
dot_8_test() ->
    ?assertEqual([[0.1925]],dot([[0.22,0.2475,0.0825]], [[0.5 ], [0.25], [0.25]])).

transpose_test() ->
    [
        ?assertEqual([[8.0]], transpose([[8.0]])),
        ?assertEqual([[1.0, 3.0], [2.0, 4.0]], transpose([[1.0, 2.0], [3.0, 4.0]])),
        ?assertEqual(
            [[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]], transpose([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        )
    ].

matmul_test() ->
    [
        ?assertEqual([[6.0]], matmul([[2.0]], [[3.0]])),
        ?assertEqual(
            [[5.0, 11.0], [11.0, 25.0]], matmul([[1.0, 2.0], [3.0, 4.0]], [[1.0, 3.0], [2.0, 4.0]])
        ),
        ?assertEqual([[30, 36, 42]], matmul([[1, 2, 3]], [[1, 2, 3], [4, 5, 6], [7, 8, 9]])),
        ?assertEqual(
            [[14], [32], [50]], matmul([[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[1], [2], [3]])
        ),
        ?_assertException(error, function_clause, matmul([[2.0], [3]], [[3.0]]))
    ].

inv_test() ->
    [
        ?assertEqual([[0.125]], inv([[8]])),
        ?assertEqual([[1.0, 0.0], [0.0, 0.5]], inv([[1, 0], [0, 2]])),
        ?assertEqual(
            [[-1.0, -1.0, 2.0], [-1.0, 0.0, 1.0], [2.0, 1.0, -2.0]],
            inv([[1, 0, 1], [0, 2, 1], [1, 1, 1]])
        )
    ].

zeros_test() ->
    [
        ?assertEqual([], zeros(0)),
        ?assertEqual([0, 0], zeros(2)),
        ?assertEqual([[]], zeros(0, 0)),
        ?assertEqual([[]], zeros(1, 0)),
        ?assertEqual([[]], zeros(0, 2)),
        ?assertEqual([[0]], zeros(1, 1)),
        ?assertEqual([[0, 0, 0], [0, 0, 0]], zeros(2, 3))
    ].

ones_test() ->
    [
        ?assertEqual([], ones(0)),
        ?assertEqual([1, 1], ones(2)),
        ?assertEqual([[]], ones(0, 0)),
        ?assertEqual([[]], ones(1, 0)),
        ?assertEqual([[]], ones(0, 2)),
        ?assertEqual([[1]], ones(1, 1)),
        ?assertEqual([[1, 1, 1], [1, 1, 1]], ones(2, 3))
    ].

fill_test() ->
    [
        ?assertEqual(zeros(2), fill(2, 0)),
        ?assertEqual(ones(3), fill(3, 1)),
        ?assertEqual([2, 2, 2, 2], fill(4, 2)),
        ?assertEqual([[]], fill(0, 0, 8)),
        ?assertEqual([[]], fill(1, 0, 7)),
        ?assertEqual([[]], fill(0, 2, 1)),
        ?assertEqual(zeros(2, 3), fill(2, 3, 0)),
        ?assertEqual(ones(3, 4), fill(3, 4, 1)),
        ?assertEqual([[2], [2], [2], [2]], fill(4, 1, 2))
    ].
