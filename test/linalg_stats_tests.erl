-module(linalg_stats_tests).
-import(linalg, [mean/1, median/1, var/1, cov/2]).
-include_lib("eunit/include/eunit.hrl").

mean_test() ->
    [
        ?assertEqual(1.0, mean([1.0])),
        ?assertEqual(1.5, mean([1.0, 2.0])),
        ?assertEqual(2.0, mean([1.0, 2.0, 3.0])),
        ?assertEqual(2.5, mean([1.0, 2.0, 3.0, 4.0])),
        ?assertEqual(12 / 5, mean([1.0, 2.0, 2.0, 3.0, 4.0]))
    ].

median_test() ->
    [
        ?assertEqual(1.0, median([1.0])),
        ?assertEqual(1.5, median([1.0, 2.0])),
        ?assertEqual(2.0, median([1.0, 2.0, 3.0])),
        ?assertEqual(2.5, median([1.0, 2.0, 3.0, 4.0])),
        ?assertEqual(2.0, median([1.0, 2.0, 2.0, 3.0, 4.0]))
    ].

var_test() ->
    [
        ?assertEqual(0, var([1.0])),
        ?assertEqual(1 / 4, var([1.0, 2.0])),
        ?assertEqual(2 / 3, var([1.0, 2.0, 3.0])),
        ?assertEqual(1.25, var([1.0, 2.0, 3.0, 4.0]))
    ].

covar_test() ->
    [
        ?assertEqual([[na, na], [na, na]], cov([1.0], [1.0])),
        ?assertEqual([[0.5, 0.5], [0.5, 0.5]], cov([1.0, 2.0], [3.0, 4.0])),
        ?assertEqual([[1.0, 1.0], [1.0, 1.0]], cov([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]))
    ].
