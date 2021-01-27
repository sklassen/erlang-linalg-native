-module(linalg_cholesky).
-vsn('1.0').
-author('simon.klassen').
-export([crout/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% example
% 25  15  -5                 5   0   0
% 15  18   0         -->     3   3   0
% -5   0  11                -1   1   3
% >>> np.linalg.cholesky([[25,15,-5],[15,18,0],[-5,0,11]])
%
% 18  22   54   42           4.24264    0.00000    0.00000    0.00000
% 22  70   86   62   -->     5.18545    6.56591    0.00000    0.00000
% 54  86  174  134          12.72792    3.04604    1.64974    0.00000
% 42  62  134  106           9.89949    1.62455    1.84971    1.39262
%>>> np.linalg.cholesky([[18, 22,  54,  42], [22, 70,  86,  62], [54, 86, 174, 134], [42, 62, 134, 106]])
%

-spec crout(matrix()) -> matrix().
crout(RowWise)->
   crout(0,RowWise).

crout(N,Matrix) when N==length(Matrix)->
     Matrix;
crout(N,Matrix)->
     % peel
     {Top,Bottom}=lists:split(N,Matrix),
     {[SL|Left],Minor}=lists:unzip([lists:split(N,Rs)|| Rs<-Bottom]),
     [[H|Row]|Rows]=Minor,
     {Col,SubMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-Rows]),
     % calc
     S=linalg:sumsq(SL),
     NewH=sqrt(H-S),
     NewRow=[0||_<-Row],
     NewCol=[[(Y-linalg:dot(SL,lists:sublist(linalg:row(J+1,[SL|Left]),N)))/NewH]||{J,[Y]}<-lists:zip(lists:seq(1,length(Col)),Col)],
     % rebuild peel
     NewBottom=[A++B++C||{A,B,C}<-lists:zip3([SL|Left],[[NewH]|NewCol],[NewRow|SubMinor])],
     crout(N+1,lists:append(Top,NewBottom)).

sqrt(X) when X<0 ->
   erlang:error({error,matrix_not_positive_defined});
sqrt(X) ->
   math:sqrt(X).


