-module(linalg_svd).
-vsn('1.0').
-author('simon.klassen').
-import(linalg,[shape/1,transpose/1,eye/1,matmul/2,norm/1,sub/2,mul/2,divide/2]).
-export([svd/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

-spec svd(matrix()) -> {matrix(),matrix(),matrix()}.
svd(RowWise)->
   case shape(RowWise) of
      {M,N} when N>M -> {error,"Need more rows than columns"};
      {M,N} -> svd(RowWise,linalg_qr:qr(RowWise),M,N,50)
   end.

svd(Matrix,U,M,N,Count)->
   {undefined,undefined,undefined}.


pytag(X,Y) ->
   case {abs(X),abs(Y)} of
           {A,B} when A>B -> A*math:sqrt(1+(B*B/A/A));
           {_,B} when B==0 ->0.0;
           {A,B} when A<B -> B*math:sqrt(1+(A*A/B/B))
   end.
