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
    {ATop,[A0|ABottom]}=lists:split(N,linalg:col(1,AMatrix)),

    {LTop,LBottom}=lists:split(N,L),
    {LLeft,[[L0|LRow]|LRows]}=lists:unzip([lists:split(N,Rs)|| Rs<-LBottom]),
    {_LCol,LMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-LRows]),

    {UTop,UBottom}=lists:split(N+1,U),
    {ULeft,[[_|URow]|URows]}=lists:unzip([lists:split(N,Rs)|| Rs<-UTop]),
    {_UCol,UMinor}=lists:unzip([lists:split(1,Rs)|| Rs<-URows]),

    LTriangle=lower(element(1,lists:split(N+1,L))),
    UNewA=u(ATop++[A0],LTriangle),

    LNewA=l([L0|ABottom],LLeft,UNewA),
    
    LNew=[A++B++C||{A,B,C}<-lists:zip3(LLeft,LNewA,[LRow|LMinor])],
    UNew=[A++B++C||{A,B,C}<-lists:zip3(ULeft,UNewA,[URow|UMinor])],

    lu(linalg:col(-1,AMatrix),{LTop++LNew,UNew++UBottom},N+1).


lower(M)->
    lower(M,0,[]).
lower([],_,Acc)->
    lists:reverse(Acc);
lower([Row|Matrix],I,Acc)->
    lower(Matrix,I+1,[element(1,lists:split(I,Row))|Acc]).
  


u([A|As],[_|Ls])->
    u(As,Ls,[A/1.0]).
u([],_Ls,Us)->
    linalg:transpose([lists:reverse(Us)]);
u([A|As],[Ls|LTail],Us)->
    S1=linalg:dot(lists:reverse(Ls),Us),
    u(As,LTail,[(A-S1)|Us]).

l([A|As],[_|Ls],Us)->
    l(As,Ls,Us,[A]).
l([],_Ls,_Us,Ls)->
    linalg:transpose([lists:reverse(Ls)]);
l([A|As],[Ls|LTail],Us,NewL)->
    UCol=element(1,lists:split(length(Ls),lists:flatten(Us))),
    [Ujj]=lists:last(Us),
    S2=linalg:dot(Ls,UCol),
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

