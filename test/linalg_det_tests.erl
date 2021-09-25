-module(linalg_det_tests).
-import(linalg,[det/1]).
-include_lib("eunit/include/eunit.hrl").

det_test()->
    [
		?assertEqual(6,det([[6]])),
		?assertEqual(-2,det([[1,2],[3,4]])),
		?assertEqual(-306,det([[6,1,1], [4,-2,5], [2,8,7]])),
		?assertEqual(2,det([[0,0,1,2], [0,0,3,4], [0,1,0,0],[1,0,0,0]])),
		?assertEqual(4,det([[1,2,0,0], [3,4,0,0], [0,0,1,2],[0,0,3,4]])),
		?assertEqual(4.0,det([[1,2,0,0,0], [3,4,0,0,0], [0,0,1,0,0], [0,0,0,1,2],[0,0,0,3,4]]))
    ].

