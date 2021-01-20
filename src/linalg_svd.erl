-module(linalg_svd).
-vsn('1.0').
-author('simon.klassen').
-import(linalg,[shape/1,transpose/1,eye/1,matmul/2,norm/1,sub/2,mul/2,divide/2]).
-import(lists,[reverse/1,merge/1,merge/2,nth/2,nthtail/2,split/2,keyreplace/4]).
-export([svd/1]).

-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

-spec svd(matrix()) -> {matrix(),matrix(),matrix()}.
svd(RowWise)->
   svd(RowWise,shape(RowWise)).

svd(_Matrix,{M,N}) when N>M ->
   {error,"Need more rows than columns"};

svd(U,{M,N})->
   Lhs = [{I,0.0} || I <- lists:seq(0,N-1)],
   Rhs = [{I,0.0} || I <- lists:seq(0,N-1)],
   svd(0,{M,N},{Lhs,Rhs},U).

svd(I,{M,N},{Lhs,Rhs},U) when I<(N-1) ->
   %io:format("Matrix (~p) ~p (G=~p)~n",[I,U,PreG]),
   {Top,Left,[[MCell|MHead]|Minor]}=minor(I,U),

   % Axis 1
   Vector1=[X||[X|_]<-[[MCell|MHead]|Minor]],
   {NewF,NewG,NewH}=fgh(Vector1),
   NewVector=[(NewF-NewG)|MHead],
   NewMinor=newaxis1(NewH,[NewVector|Minor]),

   % Axis 2
   [[Scalar2|Vector2]|Minor2]=NewMinor,
   {NextF,NextG,NextH}=fgh(Vector2),
   [_|VTail]=Vector2,
   NextVector=[Scalar2]++[NextF-NextG]++VTail,
   [_|NextMinor]=newaxis2(NextH,[NextVector|Minor2]),

   NewLhs=keyreplace(I,1,Lhs,{I,NewG}),
   NewRhs=keyreplace(I+1,1,Rhs,{I+1,NextG}),
   NewU=major({Top,Left,[NextVector|NextMinor]}),
   
   svd(I+1,{M,N},{NewLhs,NewRhs},NewU);

svd(I,{M,N},{Lhs,Rhs},U) when I==(N-1) ->
   io:format("U=~p~n",[U]),
   io:format("L=~p~n",[Lhs]),
   io:format("R=~p~n",[Rhs]),
   Left=[L||{_,L}<-Lhs],
   Right=[R||{_,R}<-Rhs],

   io:format("Right=~p~n",[right(Right,{M,N},U)]),
   {U,[L||{_,L}<-Lhs],[R||{_,R}<-Rhs]}.

right(Eigens,{M,N},U)->
   right(reverse(Eigens),{M,N},U,[[1]]).
right([_],{_,_},_,Vt)->
   Vt;
right([0.0|Tail],{M,N},U,Vt)->
   right(Tail,{M,N},U,pad(Vt));
right([E|Tail],{M,N},U,Vt)->
   I=length(Tail),
   {Top,Bottom,Minor}=minor(I,U),
   io:format("Minor=~pth ~p ~p ~p~n",[I,Top,Bottom,Minor]),
   io:format("U=~pth ~p~n",[I,U]),
   H=E*linalg:cell(I,I+1,U),
   io:format("H=~p x ~p = ~p~n",[I,linalg:cell(I,I+1,U),H]),
   right(Tail,{M,N},U,pad(Vt)).

pad(M)->
   Top=[1|[0||_<-lists:seq(1,length(M))]],
   Bottom=[[0|Row]||Row<-M],
   [Top|Bottom].

fgh([])->{0,0,0};
fgh([F|List])->
   %io:format("fgh ~p~n",[[F|List]]),
   S=linalg:sumsq([F|List]),
   G=case F<0 of true->1; false->-1 end * math:sqrt(S),
   H=F*G-S,
   {F,G,H}.

minor(N,Matrix)->
   {Top,Bottom}=split(N,Matrix),
   S=[split(N,Row)|| Row<-Bottom],
   {Left,Minor}=lists:unzip(S),
   {Top,Left,Minor}.

major({[],_,Minor})->
   Minor;
major({Top,Left,Minor})->
   Bottom=[A++B||{A,B}<-lists:zip(Left,Minor)],
   Top++Bottom.

newaxis1(_,[[X]]) ->
   [[X]];
newaxis1(H,Matrix) ->
   linalg:transpose(scaleaxis(H,linalg:transpose(Matrix))).

newaxis2(_,[[X]]) ->
   [[X]];
newaxis2(H,Matrix) ->
   {[[Scalar|Vector]],Left,Minor}=minor(1,Matrix),
   [_|NewMinor]=scaleaxis(H,[Vector|Minor]),
   major({[[Scalar|Vector]],Left,NewMinor}).

scaleaxis(H,[Vec|Matrix]) ->
   [Vec|scaleaxis(Matrix,H,Vec,[])].
scaleaxis([Col|Tail],H,Vec,Acc) ->
   %io:format("scaleaxis ~p ~p = ~p~n",[Col,Vec,lists:zip(Col,Vec)]),
   F=linalg:dot(Col,Vec)/H,
   %io:format("S ~p~n",[S]),
   %io:format("H ~p~n",[H]),
   %io:format("F ~p~n",[F]),
   NewCol=[A+F*B||{A,B}<-lists:zip(Col,Vec)],
   scaleaxis(Tail,H,Vec,[NewCol|Acc]);
scaleaxis([],_,_,Acc) ->
   reverse(Acc).


%pythag(_,0) ->
%       0;
%pythag(X,Y) ->
%   case {abs(X),abs(Y)} of
%           {A,B} when A>B -> A*math:sqrt(1+(B*B/A/A));
%           {A,B} when A=<B -> B*math:sqrt(1+(A*A/B/B))
%   end.
