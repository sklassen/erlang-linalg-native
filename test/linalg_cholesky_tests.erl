-module(linalg_cholesky_tests). 
-import(linalg_cholesky,[cholesky/1]).
-include_lib("eunit/include/eunit.hrl").

cholesky_1x1_test() ->
    ?assertEqual([[5.0]],cholesky([[25]])).

cholesky_2x2_test() ->
    ?assertEqual([[5.0,0],[3.0,3.0]],cholesky([[25,15],[15,18]])).

cholesky_3x3_test() ->
    ?assertEqual([[5.0,0,0],[3.0,3.0,0],[-1.0,1.0,3.0]],cholesky([[25,15,-5],[15,18,0],[-5,0,11]])).

cholesky_4x4_test() ->
    ?assertEqual([[4.242640687119285,0,0,0], [5.185449728701349,6.565905201197403,0,0], [12.727922061357857,3.0460384954008553,1.6497422479090704,0], [9.899494936611665,1.6245538642137891,1.849711005231382, 1.3926212476455924]],cholesky([[18, 22,  54,  42], [22, 70,  86,  62], [54, 86, 174, 134], [42, 62, 134, 106]])).
