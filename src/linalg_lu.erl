-module(linalg_lu).
-vsn('1.0.1').
-author('simon.klassen').
-export([lu/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% linalg_lu:lu([ [1.0, 3.0, 5.0], [2.0, 4.0, 7.0], [1.0, 1.0, 0.0]]).
% linalg_lu:lu([[11,9,24,2],[1,5,2,6],[3,17,18,1],[2,5,7,1]]).

pivot(RowWise)->
    pivot({linalg:eye(length(RowWise)),linalg:abs(RowWise)},[],1).

pivot({[],[]},Pivot,_I)->
    lists:reverse(Pivot);

pivot({OldPs,OldRows},NewPs,I) ->
    ArgMax=linalg:argmax(linalg:col(I,OldRows)),
    NextPs=[linalg:row(ArgMax,OldPs)|NewPs],
    pivot({linalg:row(-ArgMax,OldPs),linalg:row(-ArgMax,OldRows)},NextPs,I+1).

-spec lu(matrix()) -> {matrix(),matrix(),matrix()}.
lu(RowWise)->
    N=length(RowWise),
    Pivot=pivot(RowWise),
    {Pivot,lu(linalg:matmul(Pivot,RowWise),{linalg:eye(N),linalg:zeros(N,N)},0)}.

lu(_,{L,U},N) when N==length(L)->
    {L,U};

lu(AMatrix,{L,U},N)->
    io:format("~p) ~p~n",[N,AMatrix]),
    io:format("L) ~p~n",[L]),
    io:format("U) ~p~n",[U]),
    {ATop,[A0|ABottom]}=lists:split(N,linalg:col(1,AMatrix)),
    io:format("a) ~p;~p;~p~n",[ATop,A0,ABottom]),

    {LTop,LBottom}=lists:split(N,L),
    {LLeft,[[L0|LRow]|LRows]}=lists:unzip([lists:split(N,Rs)|| Rs<-LBottom]),
    {LCol,LMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-LRows]),

    {UBottom,UTop}=lists:split(N+1,U),
    {ULeft,[[U0|URow]|URows]}=lists:unzip([lists:split(N,Rs)|| Rs<-UBottom]),
    {UCol,UMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-URows]),

    io:format(" l) ~p;~p;~p,~p (#~p)~n",[LLeft,[L0],LRow,LCol,length(LTop)]),
    io:format(" u) ~p;~p;~p,~p (#~p)~n",[ULeft,[U0],URow,UCol,length(UBottom)]),

    LNewA=[[A0]|linalg:transpose([ABottom])],
    UNewA=linalg:transpose([ATop])++[[A0]],

    io:format("*l) ~p~n",[LNewA]),
    io:format("*u) ~p~n",[UNewA]),
    
    io:format("zl) ~p~n",[{LLeft,LNewA,[LRow|LMinor]}]),
    io:format("zu) ~p~n",[{ULeft,UNewA,[URow|UMinor]}]),
    LNew=[A++B++C||{A,B,C}<-lists:zip3(LLeft,LNewA,[LRow|LMinor])],
    UNew=[A++B++C||{A,B,C}<-lists:zip3(ULeft,UNewA,[URow|UMinor])],

    lu(linalg:col(-1,AMatrix),{LTop++LNew,UNew++UTop},N+1).

