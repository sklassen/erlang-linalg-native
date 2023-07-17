-module(linalg_roots_tests).
-import(linalg, [roots/1]).
-import(linalg_complex, [add/1, add/2, mltp/1, mltp/2, reciprocal/1]).
-define(NODEBUG, true). % Define NODEBUG for a quiet test.
-include_lib("eunit/include/eunit.hrl").
-define(SMALL, 1.0e-9).

roots_0_test() ->
    ?assert(roots([]) == []).

roots_1_test() ->
    ?assert(roots([1]) == []).

roots_01_test() ->
    ?assert(roots([0,1]) == []).

roots_2_test() ->
    ?assert(roots([1, 2]) == [-2.0]),
    Range = lists:seq(-10,10),
    ParamsList = [[A,B]||A <- Range, B <- Range, A=/=0],
    assert_roots(ParamsList).

roots_3_test() ->
    ?assert(roots([1, -3, 2]) == [1, 2]),
    Range = lists:seq(-10,10),
    ParamsList = [[A,B,C]||A <- Range, B <- Range, C <- Range],
    assert_roots(ParamsList).

roots_4_1_test() ->
    ?assert(roots([1, -6, 11, -6]) == [1, 2, 3]),
    Range = lists:seq(-10,10),
    ParamsList = [[1,B,C,D]||B <- Range, C <- Range, D <- Range],
    assert_roots(ParamsList).

roots_4_n_test() ->
    ?assertEqual([-3.0, -0.5, 2.0],[round(X*10)/10||X<-roots([2, 3, -11, -6])]),
    Range = lists:seq(-10,10),
    ParamsList = [[A,B,C,D]||A <- Range, B <- Range, C <- Range, D <- Range],
    assert_roots(ParamsList).


%% @doc Assert roots with any given params.

-type vector() :: [number()].
%spec roots(Params::vector()) -> Roots::vector().

-spec assert_roots([Params::vector()]) -> ok.
assert_roots([]) ->
    ok;
assert_roots([Params|OtherParams]) ->
    Roots = roots(Params),
    assert_roots(Params, Roots),
    assert_roots(OtherParams).

-spec assert_roots(Params::vector(), Roots::vector()) -> ok.
assert_roots(_Params, []) ->
    ok;
assert_roots(Params, [X|OtherRoots]) ->
    ?debugFmt("Params=~p, X=~p", [Params, X]),
    FoldFun = fun(Param, {XpN, Ans}) ->
        {mltp(XpN, X), add(Ans, mltp(Param, XpN))}
    end,
    Res = case lists:foldr(FoldFun, {1, 0}, Params) of
        {_, {R, I}} -> math:sqrt(R*R+I*I);
        {_, R} -> R
    end,
    ?debugFmt("abs=~p~n", [Res]),
    ?assert(abs(Res) < ?SMALL),
    assert_roots(Params, OtherRoots).
