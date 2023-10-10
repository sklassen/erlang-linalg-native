-module(linalg_poly_test).
-import(linalg_poly, [polyfit/3,polyval/2]).
-include_lib("eunit/include/eunit.hrl").
-define(PRECISION,1000000).

fit_0_test() ->
    Xs=[-1,0.0,1],
    Ys=[-25,0,25],
    ?assertEqual([0.0],approx(polyfit(Xs,Ys,0))).

fit_1_test() ->
    Xs=[-1,0.0,1],
    Ys=[55,0,25],
    ?assertEqual([26.666667,-15.0],approx(polyfit(Xs,Ys,1))).

interp_1_test() ->
    Xs=[-1,0.0,1],
    Ys=[55,0,25],
    CoEff = polyfit(Xs,Ys,1),
    ?assertEqual([41.516667,40.166667,34.166667,26.666667,19.166667,11.666667],approx([polyval(CoEff,X) || X <- [-0.99,-0.9,-0.5,0.0,0.5,1.0]])).

fit_2_test() ->
    Xs=[-1,0.0,1],
    Ys=[55,0,25],
    ?assertEqual([0.0,-15.0,40.0],approx(polyfit(Xs,Ys,2))).

interp_2_test() ->
    Xs=[-1,0.0,1],
    Ys=[55,0,25],
    CoEff = polyfit(Xs,Ys,2),
    ?assertEqual([54.054,45.9,17.5,0.0, 2.5,25.0],approx([polyval(CoEff,X) || X <- [-0.99,-0.9,-0.5,0.0,0.5,1.0]])).

approx(Fs)->
    approx(Fs,[]).
approx([],Acc)->
    lists:reverse(Acc);
approx([F|Tail],Acc)->
    approx(Tail,[round(F*?PRECISION)/?PRECISION|Acc]).
