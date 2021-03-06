-module(linalg_qr_tests). 
-import(linalg,[qr/1]).
-include_lib("eunit/include/eunit.hrl").

qr_1_test() ->
	?assertEqual({[[1]],[[-49.0]]},qr([[-49.0]])).

qr_2_test() ->
	?assertEqual({[[-0.28000000000000025,0.96],[0.96,0.28]],[[175.0,-70.0],[7.105427357601002e-15,-35.0]]},qr([[-49,-14],[168,-77]])).

qr_3_test() ->
	?assertEqual({[[0.8944271909999157,0.44721359549995804],[0.44721359549995804,-0.8944271909999155]],[[2.2360679774997894,4.919349550499537],[5.551115123125783e-16,-0.8944271909999142]]},qr([[2,4],[1,3]])).
