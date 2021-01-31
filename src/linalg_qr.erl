-module(linalg_qr).
-vsn('1.0').
-author('simon.klassen').
-export([qr/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% wikipedia
% linalg_qr:householder([[12,-51,4],[6,167,-68],[-4,24,-41]]).
% np.linalg.qr([[12,-51,4],[6,167,-68],[-4,24,-41]])
% (array([[-0.85714286,  0.39428571,  0.33142857],
%      [-0.42857143, -0.90285714, -0.03428571],
%      [ 0.28571429, -0.17142857,  0.94285714]]), 
% array([[ -14.,  -21.,   14.],
%      [   0., -175.,   70.],
%      [   0.,    0.,  -35.]]))


-spec householder(matrix()) -> {matrix(),matrix()}.
householder(ColumnWise)->
    householder(ColumnWise,{linalg:eye(length(ColumnWise)),ColumnWise}).

-spec householder(matrix(),{matrix(),matrix()}) -> {matrix(),matrix()}.
householder([[_]],{QMatrix,RMatrix})->
    {QMatrix,linalg:transpose(RMatrix)};

householder(Matrix,{QMatrix,RMatrix})->
    HMatrix=hmatrix(Matrix),
    NewQ=pad(HMatrix,length(QMatrix)-length(HMatrix)),
    Minor=minor(linalg:matmul(Matrix,HMatrix)),
    householder(Minor,{linalg:matmul(QMatrix,NewQ),linalg:matmul(RMatrix,NewQ)}).

-spec hmatrix(matrix()) -> matrix().
hmatrix([[H|Col]|_])->
    Alpha = linalg:norm([H|Col]),
    U =[(H-Alpha)|Col],
    V =linalg:divide([U],linalg:norm(U)),
    linalg:sub(linalg:eye(length(U)),linalg:mul(2,linalg:matmul(linalg:transpose(V),V))).

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
   householder(linalg:transpose(RowWise)).

