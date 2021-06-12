-module(linalg_cholesky).
-vsn('1.0').
-author('simon.klassen').
-export([cholesky/1]).

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

crout(N,Matrix) when N==length(Matrix)->
     Matrix;
crout(N,Matrix)->
     {Top,Bottom}=lists:split(N,Matrix),
     {[SL|Left],Minor}=lists:unzip([lists:split(N,Rs)|| Rs<-Bottom]),
     [[H|Row]|Rows]=Minor,
     {Col,SubMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-Rows]),
     % print
     io:format("~nGO M=~p (~p)~n",[Matrix,N]),
     io:format("H=~p~n",[H]),
     io:format("Col=~p~n",[Col]),
     io:format("Row=~p~n",[Row]),
     io:format("Left=~p~n",[Left]),
     io:format("SL=~p~n",[SL]),
     S=linalg:sumsq(SL),
     NewH=sqrt(H-S),
     NewRow=[0||_<-Row],
     NewCol=[[(Y-dot(SL,lists:sublist(linalg:row(J+1,[SL|Left]),N)))/NewH]||{J,[Y]}<-lists:zip(lists:seq(1,length(Col)),Col)],
     io:format("new col ~p ~n",[NewCol]),
     % rebuild
     NewBottom=[A++B++C||{A,B,C}<-lists:zip3([SL|Left],[[NewH]|NewCol],[NewRow|SubMinor])],
     crout(N+1,lists:append(Top,NewBottom)).

dot(A,B)->
   io:format("dot(~p,~p)=~p~n",[A,B,linalg:dot(A,B)]),
   linalg:dot(A,B).

sqrt(X) when X<0 ->
   erlang:error({error,matrix_not_positive_defined});
sqrt(X) ->
   math:sqrt(X).

-spec cholesky(matrix()) -> matrix().
cholesky(RowWise)->
   crout(0,RowWise).

