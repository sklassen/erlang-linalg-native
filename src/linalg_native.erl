-module(linalg_native). 
-vsn('1.0').
-author('simon.klassen').
-import(lists,[reverse/1,append/2,nth/2,seq/2,split/2,zip/2,foldl/3]).
-export([row/2,col/2]). 
-export([transpose/1,det/1,inv/1,shape/1,dot/2,matmul/2]). 
-export([zeros/1,ones/1,identity/1,diag/1]).
-export([exp/1,log/1]).
-export([add/2,sub/2,mul/2,divide/2,pow/2]).

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
zeros(N) -> 
	[ 0.0 ||_<-seq(1,N)].

ones(N) -> 
	[ 1.0 ||_<-seq(1,N)].

% generation (matrix)
%sequential(NR,NC) ->
%	[ [ (((R-1)*NC)+C)/1.0 || C<-seq(1,NC)] || R<-seq(1,NR)].

diag([H|_]=X) when is_number(H)->
	[ [ case R of C -> nth(R,X); _->0.0 end||R<-seq(1,length(X))] || C<-seq(1,length(X))].

identity(N) ->
	diag(ones(N)).

% Transformation
transpose([[]]) -> [];
transpose([[X]]) -> [[X]];
transpose([[] | XXs]) -> transpose(XXs);
transpose([[X | Xs] | XXs]) -> [[X | [H || [H | _Tail ] <- XXs]] | transpose([Xs | [Tail || [_|Tail] <- XXs]])].


% Sum Product
dot([],[],_) ->[];
dot([A],[B],Sum) ->Sum+A*B;
dot([A|VecA],[B|VecB],Sum) ->dot(VecA,VecB,Sum+A*B).
dot(VecA,VecB) ->dot(VecA,VecB,0).

% Matrix Multiplication
matmul(M1, M2) -> 
	Inner = length(M2),
	NCols = length(nth(1,M2)), 
	NRows = length(M1), 
	matmul(Inner, NCols, NRows,[], M1, transpose(M2)).

matmul(_, _, 0, MM, _, _) -> MM;
matmul(I, C, R, MM, M1, M2) ->
	NewRow = rowmult(I, C, R, [], M1, M2),
	matmul(I, C, R-1, [NewRow|MM], M1, M2).


rowmult(_, 0, _, L, _, _) -> L;
rowmult(I, C, R, L, M1, M2) -> 
	SumProd = dot(nth(R,M1),nth(C,M2)),
	rowmult(I, C-1, R, [SumProd|L], M1, M2).


row(I,Matrix) when I>0 ->
    [nth(I,Matrix)];
row(I,Matrix) when I<0 ->
    {A,[_|B]}=split(-(I+1),Matrix),
    append(A,B).

col(J,Matrix) ->
    transpose(row(J,transpose(Matrix))).

% Solves
det([[X]])->
    X;
det([[A,B],[C,D]])->
    A*D-B*C;
det([H|Tail])->
    foldl(fun(A,Sum)->Sum+A end,0,[pow(-1,J-1)*X*det(col(-J,Tail))||{J,X}<-zip(seq(1,length(H)),H)]).

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


% Arithmetric
exp(M)->
    sig1(M,fun(X)->math:exp(X) end,[]).

log(M)->
    sig1(M,fun(X)->math:exp(X) end,[]).

add(M1,M2)->
    sig2(M1,M2,fun(A,B)->A+B end,[]).

sub(M1,M2)->
    sig2(M1,M2,fun(A,B)->A-B end,[]).

mul(M1,M2)->
    sig2(M1,M2,fun(A,B)->A*B end,[]).

divide(M1,M2)->
    sig2(M1,M2,fun(A,B)->case B of 0-> na; B-> A/B end end,[]).

pow(M1,M2)->
    sig2(M1,M2,fun(A,B)->math:pow(A,B) end,[]).

% private functions

sig1(X,Fun,_) when is_number(X)->
    Fun(X);
sig1([],_Fun,Acc)->
    reverse(Acc);
sig1([H|_]=Vector,Fun,[]) when is_number(H)->
    [Fun(X)||X <-Vector];
sig1([R1|Matrix],Fun,Acc)->
    sig1(Matrix,Fun,[[Fun(X)||X <-R1]|Acc]).


sig2([],[],_Fun,Acc)->
    reverse(Acc);
sig2(X,[],_Fun,Acc) when is_number(X)->
    reverse(Acc);
sig2([],X,_Fun,Acc) when is_number(X)->
    reverse(Acc);
sig2(A,B,Fun,[]) when is_number(A) andalso is_number(B)->
    Fun(A,B);
sig2(A,[B|Vector],Fun,[]) when is_number(A) andalso is_number(B)->
    [Fun(A,X)||X<-[B|Vector]];
sig2([A|Vector],B,Fun,[]) when is_number(A) andalso is_number(B)->
    [Fun(X,B)||X<-[A|Vector]];
sig2([A|VectorA],[B|VectorB],Fun,[]) when is_number(A) andalso is_number(B)->
    [Fun(X,Y)||{X,Y}<-zip([A|VectorA],[B|VectorB])];
sig2(A,[R2|M2],Fun,Acc) when is_number(A)->
    sig2(A,M2,Fun,[[Fun(A,B)||B<-R2]|Acc]);
sig2([R1|M1],B,Fun,Acc) when is_number(B)->
    sig2(M1,B,Fun,[[Fun(A,B)||A<-R1]|Acc]);
sig2([R1|M1],[R2|M2],Fun,Acc)->
    sig2(M1,M2,Fun,[[Fun(A,B)||{A,B}<-zip(R1,R2)]|Acc]).
