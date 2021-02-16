-module(linalg_lu).
-vsn('1.0.1').
-author('simon.klassen').
-export([swap/3,lu/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% linalg_lu:lu([[1.0,3.0,5.0],[2.0,4.0,7.0],[1.0,1.0,0.0]]).
% linalg_lu:lu([[11,9,24,2],[1,5,2,6],[3,17,18,1],[2,5,7,1]]).


-spec lu(matrix()) -> {matrix(),matrix(),matrix()}.
lu(RowWise)->
    N=length(RowWise),
    Pivot=pivot(RowWise),
    {L,U}=lu(linalg:matmul(Pivot,RowWise),{linalg:eye(N),linalg:zeros(N,N)},0),
    {Pivot,L,U}.

lu(_,{L,U},N) when N==length(L)->
    {L,U};

lu(AMatrix,{L,U},N)->
    io:format("~p) ~p~n",[N,AMatrix]),
    io:format("L) ~p~n",[L]),
    io:format("U) ~p~n",[U]),
    {ATop,[A0|ABottom]}=lists:split(N,linalg:col(1,AMatrix)),
    %io:format("a) ~p;~p;~p~n",[ATop,A0,ABottom]),

    {LTop,LBottom}=lists:split(N,L),
    {LLeft,[[L0|LRow]|LRows]}=lists:unzip([lists:split(N,Rs)|| Rs<-LBottom]),
    {_LCol,LMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-LRows]),

    {UTop,UBottom}=lists:split(N+1,U),
    {ULeft,[[_|URow]|URows]}=lists:unzip([lists:split(N,Rs)|| Rs<-UTop]),
    {_UCol,UMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-URows]),

    %io:format(" l) ~p;~p;~p,~p (~p//~p)~n",[LLeft,[L0],LRow,LCol,LTop,LBottom]),
    %io:format(" u) ~p;~p;~p,~p (~p/~p)~n",[ULeft,[U0],URow,UCol,UTop,UBottom]),
    
    LTriangle=lower(element(1,lists:split(N+1,L))),
    io:format("ltri) ~p~n",[LTriangle]),

    %UNewA=linalg:transpose([As]),
    UNewA=u(ATop++[A0],LTriangle),
    io:format("U*) ~p~n",[UNewA]),

    %LNewA=[[A0]|linalg:transpose([ABottom])],
    LNewA=l([L0|ABottom],LLeft,UNewA),

    %io:format("*l) ~p~n",[LNewA]),
    %io:format("*u) ~p~n",[UNewA]),
    
    %io:format("zl) ~p~n",[{LLeft,LNewA,[LRow|LMinor]}]),
    %io:format("zu) ~p~n",[{ULeft,UNewA,[URow|UMinor]}]),
    LNew=[A++B++C||{A,B,C}<-lists:zip3(LLeft,LNewA,[LRow|LMinor])],
    UNew=[A++B++C||{A,B,C}<-lists:zip3(ULeft,UNewA,[URow|UMinor])],

    lu(linalg:col(-1,AMatrix),{LTop++LNew,UNew++UBottom},N+1).

% Ragged lower triangle (without diag) tril(M,k=-1)
lower(M)->
    lower(M,0,[]).
lower([],_,Acc)->
    lists:reverse(Acc);
lower([Row|Matrix],I,Acc)->
    lower(Matrix,I+1,[element(1,lists:split(I,Row))|Acc]).
  


u([A|As],[_|Ls])->
    io:format("Ustart> A=~p L=~p~n",[[A|As],Ls]),
    u(As,Ls,[A/1.0]).
u([],_Ls,Us)->
    io:format("Uend> ~p~n",[lists:reverse(Us)]),
    linalg:transpose([lists:reverse(Us)]);
u([A|As],[Ls|LTail],Us)->
    %io:format("> A ~p rows ~p cols ~p~n",[A,Ls,Us]),
    S1=linalg:dot(lists:reverse(Ls),Us),
    %io:format("Uij ~p = (A-S1) = (~p-~p) (~p)~n",[A-S1,A,S1,lists:zip(lists:reverse(Ls),Us)]),
    u(As,LTail,[(A-S1)|Us]).

l([A|As],[_|Ls],Us)->
    io:format("Lstart> A=~p L=~p U=~p~n",[[A|As],Ls,Us]),
    l(As,Ls,Us,[A]).
l([],_Ls,_Us,Ls)->
    io:format("Lend> ~p ~n",[lists:reverse(Ls)]),
    linalg:transpose([lists:reverse(Ls)]);
l([A|As],[Ls|LTail],Us,NewL)->
    io:format("Have > A=~p L=~p U=~p ~n",[A,[Ls|LTail],Us]),
    UCol=element(1,lists:split(length(Ls),lists:flatten(Us))),
    io:format("Have > A=~p L=~p U=~p UCol=~p~n",[A,[Ls|LTail],Us,UCol]),
    [Ujj]=lists:last(Us),
    %io:format("> A=~p L=~p U=~p Ujj=~p~n",[A,Ls,UCol,Ujj]),
    %{Col,_}=lists:split(length(NewL),Us),
    %{Row,_}=lists:split(length(NewL),Ls),
    %io:format("> rows ~p cols ~p~n",[[Ls|LTail],Us]),
    S2=linalg:dot(Ls,UCol),
    io:format("L ~p = (A-S2)/Ujj = (~p-~p)/~p (~p)~n",[(A-S2)/Ujj,A,S2,Ujj,lists:zip(lists:reverse(Ls),UCol)]),
    l(As,LTail,Us,[(A-S2)/Ujj|NewL]).

swap(RowWise,I,I) ->
    RowWise;
swap(RowWise,I,J) when J<I->
    swap(RowWise,J,I);
swap(RowWise,I,J) ->
    {A,[C|Tail]}=lists:split(I-1,RowWise),
    {B,D}=lists:split(J-I,Tail),
    lists:append([A,B,[C],D]).

pivot(RowWise)->
    pivot(linalg:abs(RowWise),linalg:eye(length(RowWise)),1).
pivot([],Pivot,_I)->
    Pivot;
pivot([R|Rows],Ps,I) ->
    pivot(Rows,swap(Ps,I,I+(linalg:argmax(linalg:col(I,[R|Rows]))-1)),I+1).

