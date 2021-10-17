-module(linalg_svd_power).
-author('simon.klassen').
-import(linalg,[t/1]).
-export([svd/1]).

-define(EPSILON, 1.0e-5).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

-spec svd(matrix()) -> {matrix(),matrix(),matrix()}.
svd(RowWise)->
  Seed=rand:seed(exsplus, {42, 123534, 345345}),
  svd(RowWise,Seed,linalg:shape(RowWise)).

svd(A,Seed,{N,N}) ->
  {Q,_} = linalg:qr(linalg:random(N,N,[{seed,Seed}])),
  {Left,Right}=power_iteration(A,Q,100),
  S = linalg:diag(Right),
  {
   t(Left),
   S,
   Left
  }.


power_iteration(_A,_Q,0)->
  {error,"no convergence"};
power_iteration(A,PrevQ,N)->
  {Q,R}=linalg:qr(linalg:matmul(A,PrevQ)),
  case linalg:sum(linalg:pow(linalg:sub(Q,PrevQ),2)) of
    Small when Small < ?EPSILON -> {Q,R};
    Large -> io:format("Loop ~p~n",[Large]),power_iteration(A,Q,N-1)
  end.


