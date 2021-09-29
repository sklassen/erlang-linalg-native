-module(linalg_roots_tests).
-import(linalg, [roots/1]).
-include_lib("eunit/include/eunit.hrl").

roots_0_test() ->
    ?assert(roots([]) == []).

roots_1_test() ->
    ?assert(roots([1]) == []).

roots_2_test() ->
    ?assert(roots([1, 2]) == [-2.0]).
