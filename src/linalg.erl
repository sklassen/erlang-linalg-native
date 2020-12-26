-module(linalg). 
-vsn('1.0').
-author('simon.klassen').
-import(lists,[append/2,nth/2,seq/2,split/2,zip/2,foldl/3]).
-import(linalg_arithmetric,[mul/2,divide/2,pow/2]).
-export([row/2,col/2,cell/3]). 
-export([transpose/1,det/1,inv/1,shape/1,dot/2,matmul/2]). 
-export([zeros/1,ones/1,sequential/1,random/1,zeros/2,ones/2,sequential/2,random/2]).
-export([identity/1,diag/1,eye/1,eye/2]).

-type dim() :: non_neg_integer().
-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% linalg shape 
shape(X) when is_number(X)->
   {};
shape([X|_]=Vector) when is_number(X)->
   {length(Vector)};
shape([[X|_]|_]=Matrix) when is_number(X)->
   NRows = length(Matrix), 
   NCols = length(nth(1,Matrix)), 
   {NRows,NCols}.

% generation (vector)
-spec zeros(dim())->vector().
zeros(0) -> 
    [[]];
zeros(N) -> 
	[ 0.0 ||_<-seq(1,N)].

-spec ones(dim())->vector().
ones(0) -> 
    [[]];
ones(N) -> 
	[ 1.0 ||_<-seq(1,N)].

-spec sequential(dim())->vector().
sequential(0) ->
    [[]];
sequential(N) ->
	[ X||X<-seq(1,N)].

-spec random(dim())->vector().
random(0) ->
    [[]];
random(N) ->
	[ rand:uniform() ||_<-seq(1,N)].

% generation (matrix)
-spec zeros(dim(),dim())->matrix().
zeros(NR,NC) ->
	[ [ 0 || _<-seq(1,NC)] || _<-seq(1,NR)].

-spec ones(dim(),dim())->matrix().
ones(NR,NC) ->
	[ [ 1 || _<-seq(1,NC)] || _<-seq(1,NR)].

-spec sequential(dim(),dim())->matrix().
sequential(NR,NC) ->
	[ [ (((R-1)*NC)+C)/1.0 || C<-seq(1,NC)] || R<-seq(1,NR)].

-spec random(dim(),dim())->matrix().
random(NR,NC) ->
	[ [ rand:uniform() || _<-seq(1,NC)] || _<-seq(1,NR)].

-spec eye(dim())->matrix().
eye(0)->
     eye([[]]);
eye(N)->
     eye(N,N).

-spec eye(dim(),dim())->matrix().
eye(N,M) ->
      [ [ case {R,C} of {C,R} -> 1.0; _->0.0 end||R<-seq(1,M)] || C<-seq(1,N)].

-spec diag(vector()|matrix())->matrix()|vector().
% (V)ector V->M
diag([X|_]=V) when is_number(X)->
	[ [ case R of C -> nth(R,V); _->0.0 end||R<-seq(1,length(V))] || C<-seq(1,length(V))];

% (M)atrix M->V
diag([[X|_]|_]=M) when is_number(X)->
	[ nth(R,nth(R,M))||R<-seq(1,length(M))].

-spec identity(dim())->matrix().
identity(N) ->
	diag(ones(N)).

% Transformation
-spec transpose(matrix())->matrix().
transpose([[]]) -> [];
transpose([[X]]) -> [[X]];
transpose([[] | XXs]) -> transpose(XXs);
transpose([[X | Xs] | XXs]) -> [[X | [H || [H | _Tail ] <- XXs]] | transpose([Xs | [Tail || [_|Tail] <- XXs]])].

% Sum Product
-spec dot(vector(),vector())->scalar().
dot([],[],Sum) ->Sum;
dot([A],[B],Sum) ->Sum+A*B;
dot([A|VecA],[B|VecB],Sum) ->dot(VecA,VecB,Sum+A*B).
dot(VecA,VecB) ->dot(VecA,VecB,0).

% Matrix Multiplication
-spec matmul(matrix(),matrix())->matrix().
matmul(M1 = [H1|_], M2) when length(H1) =:= length(M2) ->
    matmul(M1, transpose(M2), []).

matmul([], _, R) -> lists:reverse(R);
matmul([Row|Rest], M2, R) ->
    matmul(Rest, M2, [rowmult(Row, M2, [])|R]).

rowmult(_, [], R) -> lists:reverse(R);
rowmult(Row, [Col|Rest], R) ->
    rowmult(Row, Rest, [dot(Row, Col)|R]).

% Note: much slower but succient. 
% matmul_zipwith(M1,M2) -> 
%   [ [ foldl(fun(X,Sum)->Sum+X end,0,lists:zipwith(fun(X,Y)->X*Y end,A,B))|| A <- transpose(M1) ]|| B <- M2 ].


% Reference
-spec row(dim(),matrix())->vector().
row(I,Matrix) when I>0 ->
    [nth(I,Matrix)];
row(I,Matrix) when I<0 ->
    {A,[_|B]}=split(-(I+1),Matrix),
    append(A,B).

-spec col(dim(),matrix())->vector().
col(J,Matrix) ->
    transpose(row(J,transpose(Matrix))).

-spec cell(dim(),dim(),matrix())->vector().
cell(I,J,Matrix) ->
    nth(J,nth(1,row(I,Matrix))).

% Solves
-spec det(matrix())->scalar().
det([[X]])->
    X;
det([[A,B],[C,D]])->
    A*D-B*C;
det([H|Tail])->
    foldl(fun(A,Sum)->Sum+A end,0,[pow(-1,J-1)*X*det(col(-J,Tail))||{J,X}<-zip(seq(1,length(H)),H)]).

-spec inv(matrix())->matrix().
inv([[X]])->
    [[1.0/X]];
inv([[A,B],[C,D]])->
    case det([[A,B],[C,D]]) of
       0.0->err;
       Det->[[D/Det,-1/Det*B],[-1/Det*C,A/Det]]
    end;
inv(M)->
    case det(M) of
       0.0->err;
       Det->divide(transpose(mul(minors(M),cofactors(M))),Det)
    end.

minors(Matrix)->
    {NRows,NCols}=shape(Matrix),
    [[det(col(-J,row(-I,Matrix)))||J<-seq(1,NCols)]||I<-seq(1,NRows)].

cofactors(Matrix)->
    {NRows,NCols}=shape(Matrix),
    [[pow(-1,I)*pow(-1,J)||J<-seq(0,NCols-1)]||I<-seq(0,NRows-1)].


