-module(linalg_reduction_tests).
-import(linalg, [sum/1, prod/1, norm/1, min/1, max/1]).
-include_lib("eunit/include/eunit.hrl").

sum_test() ->
    [
        ?assertEqual(1, sum(1)),
        ?assertEqual(6, sum([1, 2, 3])),
        ?assertEqual(12, sum([[1, 2, 3], [1, 2, 3]]))
    ].

prod_test() ->
    [
        ?assertEqual(1, prod(1)),
        ?assertEqual(6, prod([1, 2, 3])),
        ?assertEqual(36, prod([[1, 2, 3], [1, 2, 3]]))
    ].

min_test() ->
    [
        ?assertEqual(1, min(1)),
        ?assertEqual(1, min([1, 2, 3])),
        ?assertEqual(1, min([[1, 2, 3], [1, 2, 3]]))
    ].

max_test() ->
    [
        ?assertEqual(1, max(1)),
        ?assertEqual(3, max([1, 2, 3])),
        ?assertEqual(3, max([[1, 2, 3], [1, 2, 3]]))
    ].

norm_test() ->
    [
        ?assertEqual(1, norm(1)),
        ?assertEqual(math:sqrt(14), norm([1, 2, 3])),
        ?assertEqual(math:sqrt(28), norm([[1, 2, 3], [1, 2, 3]]))
    ].
