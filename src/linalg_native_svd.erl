-module(linalg_native_svd).
-vsn('1.0').
-author('simon.klassen').
-import(lists,[reverse/1,append/2,nth/2,seq/2,split/2,zip/2,foldl/3]).
-import(linalg_native,[row/2,col/2,cell/3]).
-import(linalg_native,[transpose/1,det/1,inv/1,shape/1,dot/2,matmul/2]).
-import(linalg_native,[zeros/1,ones/1,identity/1,diag/1,eye/1,eye/2]).
-import(linalg_arithmetric,[sum/1]).
-import(linalg_arithmetric,[exp/1,log/1]).
-import(linalg_arithmetric,[add/2,sub/2,mul/2,divide/2,pow/2]).
-export([h/1,householder/2,qr/1,pad/2,norm/1]).
-define(EPSILON,1/1000000).

% wikipedia
% linalg_native_svd:householder([[12,-51,4],[6,167,-68],[-4,24,-41]]).
% np.linalg.qr([[12,-51,4],[6,167,-68],[-4,24,-41]])
% (array([[-0.85714286,  0.39428571,  0.33142857],
%      [-0.42857143, -0.90285714, -0.03428571],
%      [ 0.28571429, -0.17142857,  0.94285714]]), 
% array([[ -14.,  -21.,   14.],
%      [   0., -175.,   70.],
%      [   0.,    0.,  -35.]]))


epsilon(X) when is_number(X) -> 
   case (abs(X)<?EPSILON) of 
      true->0;
      false->X 
   end. 

norm([H|_]=Vector) when is_number(H)->
    math:sqrt(sum(pow(Vector,2))).

sign(Zero) when Zero==0 orelse Zero==0.0 -> 0;
sign(Number) -> Number/abs(Number).

unit(N) ->
    [ [ case {R,C} of {1,1} -> 1.0; _->0.0 end||R<-seq(1,N)] || C<-seq(1,N)].

pad(ColWise,0)->
    ColWise;
pad(ColWise,N)->
    Row=[1|[0||_<-seq(1,length(nth(1,ColWise)))]],
    pad([Row|[[0|X]||X<-ColWise]],N-1).

minor([_|Matrix])->
    [Tail||[_|Tail]<-Matrix].


h([[H|Col]|_])->
    Alpha = norm([H|Col]),
    U =[(H-Alpha)|Col],
    V =divide([U],norm(U)),
    sub(eye(length(U)),mul(2,matmul(transpose(V),V))).

householder(ColumnWise)->
    householder(ColumnWise,{eye(length(ColumnWise)),ColumnWise}).

householder([[_]],{QMatrix,RMatrix})->
    {QMatrix,transpose(RMatrix)};

householder(Matrix,{QMatrix,RMatrix})->
    %N=1+length(QMatrix)-length(Matrix),
    %io:format("A~p:~p~n",[N,Matrix]),
    HMatrix=h(Matrix),
    %io:format("H~p:~p~n",[N,HMatrix]),
    Padded=pad(HMatrix,length(QMatrix)-length(HMatrix)),
    %io:format("P~p:~p~n",[N,Padded]),
    Minor=minor(matmul(Matrix,HMatrix)),
    %io:format("Q~p:~p~n",[N,matmul(QMatrix,Padded)]),
    %io:format("R~p:~p~n",[N,matmul(RMatrix,Padded)]),
    householder(Minor,{matmul(QMatrix,Padded),matmul(RMatrix,Padded)}).

qr(RowWise)->
   householder(transpose(RowWise)).

