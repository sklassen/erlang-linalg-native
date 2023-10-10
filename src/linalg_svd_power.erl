-module(linalg_svd_power).
-author('simon.klassen').
-import(linalg,[t/1]).
-export([svd/1]).

-define(EPSILON, 1.0e-12).
-define(SEED,rand:seed(exsplus, {42, 123534, 345345})).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

-spec svd(matrix()) -> {matrix(),matrix(),matrix()}.
svd(RowWise)->
  svd(RowWise,linalg:shape(RowWise)).

svd(A,{N,N}) ->
  {Q,_} = linalg:qr(linalg:random(N,N,[{seed,?SEED}])),
  {Left,Right}=power_iteration(A,Q,100),
  S = linalg:diag(Right),
  {
   t(Left),
   S,
   Left
  };

svd(A,{N,M}) when N < M->
  {Q,_} = linalg:qr(linalg:random(N,N,[{seed,?SEED}])),
  {Left,Right}=power_iteration(linalg:matmul(A,t(A)),Q,100),
  S=linalg:sqrt(linalg:diag(Right)),
  {
    t(Left),
    S,
    linalg:matmul(linalg:inv(linalg:diag(S)),linalg:matmul(t(Left),A))
  };

svd(A,{N,M}) when N > M->
  {Q,_} = linalg:qr(linalg:random(M,M,[{seed,?SEED}])),
  {Left,Right}=power_iteration(linalg:matmul(t(A),A),Q,100),
  S=linalg:sqrt(linalg:diag(Right)),
  {
    linalg:matmul(linalg:matmul(A,t(Left)),linalg:inv(linalg:diag(S))),
    S,
    t(Left)
  }.


power_iteration(_A,_Q,0)->
  erlang:error("linalg:svd no convergence");
power_iteration(A,PrevQ,N)->
  {Q,R}=linalg:qr(linalg:matmul(A,PrevQ)),
  case linalg:sum(linalg:pow(linalg:sub(Q,PrevQ),2)) of
    Small when Small < ?EPSILON -> {Q,R};
    _Large -> power_iteration(A,Q,N-1)
  end.


