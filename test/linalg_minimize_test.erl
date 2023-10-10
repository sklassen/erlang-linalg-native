-module(linalg_minimize_test).
-import(linalg_minimize, [minimize/3]).
-define(F,fun(X) -> X*X*X+X*X-3*X-3 end).
-define(MIN,0.721).
-include_lib("eunit/include/eunit.hrl").

minimizie_golden_test() ->
    ?assertEqual(?MIN,round(minimize(?F,{-2,2},golden)*1000)/1000).

minimizie_brute_test() ->
    ?assertEqual(?MIN,round(minimize(?F,{0.0,1.0},brute)*1000)/1000).

