-module(linalg_complex_tests).
-import(linalg_complex, [sqrt/1, qbrt/1]).
-import(linalg_complex, [add/1, add/2, mltp/1, mltp/2, reciprocal/1]).
-define(NODEBUG, true). % Define NODEBUG for a quiet test.
-include_lib("eunit/include/eunit.hrl").
-define(PRECISION,1000000).

complex_sqrt_1_test() ->
    ?assertEqual(0, sqrt(0)),
    ?assertEqual(0, sqrt({0,0})),
    ?assertEqual(2.0, sqrt(4)),
    ?assertEqual(2.0, sqrt({4,0})),
    ?assertEqual({0,2.0}, sqrt(-4)),
    ?assertEqual({0,2.0}, sqrt({-4,0})).

complex_qbrt_1_test() ->
    ?assertEqual(0, qbrt(0)),
    ?assertEqual(0, qbrt({0,0})),
    ?assertEqual(2.0, qbrt(8)),
    ?assertEqual(2.0, qbrt({8,0})),
    ?assertEqual([
        -1.709976, % -5
        {1.480883, 0.854988}, % {0,5}
        {1.480883, -0.854988}, % {0,-5}
        {1.451857, 0.493404}, % {2,3}
        {1.153228, 1.010643}, % {-2,3}
        {1.451857, -0.493404}, % {2,-3}
        {1.153228, -1.010643} % {-2,-3}
    ], approx([
        qbrt(-5),
        qbrt({0,5}),
        qbrt({0,-5}),
        qbrt({2,3}),
        qbrt({-2,3}),
        qbrt({2,-3}),
        qbrt({-2,-3})
    ])).

complex_add_1_test() ->
    ?assertEqual({7,8}, add([1, {2,0}, {0,3}, {4,5}])).

complex_add_2_test() ->
    ?assertEqual(3, add(1, 2)),
    ?assertEqual(3, add({1,0}, 2)),
    ?assertEqual(3, add(1, {2,0})),
    ?assertEqual(3, add({1,0}, {2,0})),
    ?assertEqual({3,9}, add({1,4}, {2,5})),
    ?assertEqual({-4,2}, add({1,4}, {-5,-2})),
    ?assertEqual({6.8,11.2}, add({1.2,3.4}, {5.6,7.8})).

complex_mltp_1_test() ->
    ?assertEqual({-196,83}, mltp([{2,3},{4,5},{6,7}])).

complex_mltp_2_test() ->
    ?assertEqual(6, mltp(3, 2)),
    ?assertEqual(6, mltp({3,0}, 2)),
    ?assertEqual(6, mltp(3, {2,0})),
    ?assertEqual(6, mltp({3,0}, {2,0})),
    ?assertEqual({-14, 23}, mltp({3,4}, {2,5})),
    ?assertEqual([{-15.55, 26.4}], approx([mltp({3.1,4.2}, {2.3,5.4})])).

complex_reciprocal_2_test() ->
    ?assertEqual(
        [0.5, 0.5, 0.4, 0.4, {0.153846, -0.230769}, {0.160932, 0.172031}],
        approx([
            reciprocal(2),
            reciprocal({2, 0}),
            reciprocal(2.5),
            reciprocal({2.5, 0}),
            reciprocal({2, 3}),
            reciprocal({2.9, -3.1})
        ])
    ).

complex_min_test() ->
    ?assertEqual({-10,-10}, linalg_complex:min([
        -1, 0, 1, {-10,-10}, {-10,0}, {-10,10}, {0,-10}, {0,0}, 
        {10,-10}, {10,0}, {10,10}, 1, 0, -1, {0,0}
    ])).

complex_max_test() ->
    ?assertEqual({10,10}, linalg_complex:max([
        -1, 0, 1, {-10,-10}, {-10,0}, {-10,10}, {0,-10}, {0,0}, 
        {10,-10}, {10,0}, {10,10}, 1, 0, -1, {0,0}
    ])).

approx(Fs) ->
    approx(Fs,[]).
approx([],Acc) ->
    lists:reverse(Acc);
approx([{Real,Imaginary}|Tail],Acc) ->
    approx(Tail,[{
        round(Real*?PRECISION)/?PRECISION,
        round(Imaginary*?PRECISION)/?PRECISION
    }|Acc]);
approx([Real|Tail],Acc) ->
    approx(Tail,[round(Real*?PRECISION)/?PRECISION|Acc]).
