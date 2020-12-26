-module(linalg_svd).
-vsn('1.0').
-author('simon.klassen').
-import(linalg,[transpose/1,eye/1,matmul/2]).
-import(linalg_arithmetric,[norm/1,sub/2,mul/2,divide/2]).
-export([qr/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% wikipedia
% linalg_native_svd:householder([[12,-51,4],[6,167,-68],[-4,24,-41]]).
% np.linalg.qr([[12,-51,4],[6,167,-68],[-4,24,-41]])
% (array([[-0.85714286,  0.39428571,  0.33142857],
%      [-0.42857143, -0.90285714, -0.03428571],
%      [ 0.28571429, -0.17142857,  0.94285714]]), 
% array([[ -14.,  -21.,   14.],
%      [   0., -175.,   70.],
%      [   0.,    0.,  -35.]]))


-spec householder(matrix()) -> {matrix(),matrix()}.
householder(ColumnWise)->
    householder(ColumnWise,{eye(length(ColumnWise)),ColumnWise}).

-spec householder(matrix(),{matrix(),matrix()}) -> {matrix(),matrix()}.
householder([[_]],{QMatrix,RMatrix})->
    {QMatrix,transpose(RMatrix)};

householder(Matrix,{QMatrix,RMatrix})->
    %N=1+length(QMatrix)-length(Matrix),
    %io:format("A~p:~p~n",[N,Matrix]),
    HMatrix=hmatrix(Matrix),
    %io:format("H~p:~p~n",[N,HMatrix]),
    NewQ=pad(HMatrix,length(QMatrix)-length(HMatrix)),
    %io:format("P~p:~p~n",[N,NewQ]),
    Minor=minor(matmul(Matrix,HMatrix)),
    %io:format("Q~p:~p~n",[N,matmul(QMatrix,NewQ)]),
    %io:format("R~p:~p~n",[N,matmul(RMatrix,NewQ)]),
    householder(Minor,{matmul(QMatrix,NewQ),matmul(RMatrix,NewQ)}).

-spec hmatrix(matrix()) -> matrix().
hmatrix([[H|Col]|_])->
    Alpha = norm([H|Col]),
    U =[(H-Alpha)|Col],
    V =divide([U],norm(U)),
    sub(eye(length(U)),mul(2,matmul(transpose(V),V))).

-spec minor(matrix()) -> matrix().
minor([_|Matrix])->
    [Tail||[_|Tail]<-Matrix].

-spec pad(matrix(),integer()) -> matrix().
pad(ColWise,0)->
    ColWise;
pad(ColWise,N)->
    Row=[1|[0||_<-lists:seq(1,length(lists:nth(1,ColWise)))]],
    pad([Row|[[0|X]||X<-ColWise]],N-1).

-spec qr(matrix()) -> {matrix(),matrix()}.
qr(RowWise)->
   householder(transpose(RowWise)).
