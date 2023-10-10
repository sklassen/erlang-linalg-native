-module(linalg_svd_tests).
-import(linalg, [matmul/2,diag/1,svd/1]).
-include_lib("eunit/include/eunit.hrl").

svd_2x2_test() ->
  A=[[3,2],[2,3]],
  {U,S,Vt}=svd(A),
  [
   ?assertEqual([[0.707,0.707],[0.707,-0.707]],linalg:around(U,3)),
   ?assertEqual([5, 1], linalg:around(S)),
   ?assertEqual([[0.707,0.707],[0.707,-0.707]],linalg:around(Vt,3)),
   ?assertEqual(A,linalg:around(matmul(matmul(U,diag(S)),Vt)))
  ].

svd_2x3_test() ->
  A=[[3,2,2],[2,3,-2]],
  {U,S,Vt}=svd(A),
  [
   ?assertEqual([[0.707,0.707],[0.707,-0.707]],linalg:around(U,3)),
   ?assertEqual([5, 3], linalg:around(S)),
   ?assertEqual([[0.707,0.707,0.0],[0.236,-0.236,0.943]],linalg:around(Vt,3)),
   ?assertEqual(A,linalg:around(matmul(matmul(U,diag(S)),Vt)))
  ].


svd_3x2_test() ->
  A=[[3,2],[2,3],[2,-2]],
  {U,S,Vt}=svd(A),
  [
   ?assertEqual([[0.707,0.236],[0.707,-0.236],[0.0,0.943]],linalg:around(U,3)),
   ?assertEqual([5, 3], linalg:around(S)),
   ?assertEqual([[0.707,0.707],[0.707,-0.707]],linalg:around(Vt,3)),
   ?assertEqual(A,linalg:around(matmul(matmul(U,diag(S)),Vt)))
  ].

