-module(linalg_svd_tests). 
-import(linalg,[svd/1]).
-include_lib("eunit/include/eunit.hrl").

svd_1_test() ->
    ?assertEqual({undefined,undefined,undefined},svd([[2,4],[1,3]])).
