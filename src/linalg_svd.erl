-module(linalg_svd).
-vsn('1.0').
-author('simon.klassen').
-import(linalg,[shape/1,transpose/1,eye/1,matmul/2,norm/1,sub/2,mul/2,divide/2]).
-import(lists,[reverse/1,merge/1,merge/2,nth/2,nthtail/2,split/2]).
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
   Lhs = [0.0 || _ <- lists:seq(1,N)],
   Rhs = [0.0 || _ <- lists:seq(1,N)],
   svd(0,{M,N},0,{Lhs,Rhs},U).

svd(I,{M,N},PreG,{Lhs,Rhs},U) when I<(N-1) ->
   io:format("Matrix (~p) ~p (G=~p)~n",[I,U,PreG]),
   {QHead,[_Q|QTail]}=split(I,Lhs),
   {EHead,[_E|ETail]}=split(I,Rhs),
   {Top,Left,[[MCell|MHead]|Minor]}=minor(I,U),
   [C|Col]=[X||[X|_]<-[[MCell|MHead]|Minor]],

   {NewF,NewG,NewH}=fgh([C|Col]),
   NewHead=[(NewF-NewG)|MHead],
   %io:format("NewMinor ~p ~p~n",[NewH,Minor]),
   NewMinor=newaxis1(NewH,[NewHead|Minor]),
   %io:format("NewMinor ~p~n",[NewMinor]),

   [[Scalar2|Vector2]|Minor2]=NewMinor,
   {NextF,NextG,NextH}=fgh(Vector2),
   [_|VTail]=Vector2,
   Vector3=[Scalar2]++[NextF-NextG]++VTail,
   %io:format("New One ~p~n",[Vector3]),

   %io:format("New Matrix~p~n",[[Vector3|Minor2]]),
   [_|NewAxis]=newaxis2(NextH,[Vector3|Minor2]),
   %io:format("New Axis ~p~n",[NewAxis]),

   NewLhs=QHead++[NewG|QTail],
   NewRhs=EHead++[PreG|ETail],
   NewU=major({Top,Left,[Vector3|NewAxis]}),
   %io:format("U=~p~n",[NewU]),
   %io:format("E=~p~n",[NewE]),
   %io:format("Q=~p~n",[NewQ]),
   svd(I+1,{M,N},NextG,{NewLhs,NewRhs},NewU);

svd(I,{_,N},G,{Lhs,Rhs},U) when I==(N-1) ->
   {EHead,[_E|ETail]}=split(I,Rhs),
   NewRhs=EHead++[G|ETail],
   {U,Lhs,NewRhs}.

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
   %io:format("newaxis inner ~p~n",[Matrix]),
   %io:format("newaxix minor0 ~p~n",[minor(1,Matrix)]),
   {[[Scalar|Vector]],Left,Minor}=minor(1,Matrix),
   %io:format("newaxis minor~p~n",[[Vector|Minor]]),
   [_|NewMinor]=scaleaxis(H,[Vector|Minor]),
   %io:format("newaxis output ~p~n",[NewMinor]),
   major({[[Scalar|Vector]],Left,NewMinor}).

scaleaxis(H,[Vec|Matrix]) ->
   [Vec|scaleaxis(Matrix,H,Vec,[])].
scaleaxis([Col|Tail],H,Vec,Acc) ->
   %io:format("scaleaxis ~p ~p = ~p~n",[Col,Vec,lists:zip(Col,Vec)]),
   S=linalg:dot(Col,Vec),
   F = S/H,
   %io:format("S ~p~n",[S]),
   %io:format("H ~p~n",[H]),
   %io:format("F ~p~n",[F]),
   NewCol=[A+F*B||{A,B}<-lists:zip(Col,Vec)],
   scaleaxis(Tail,H,Vec,[NewCol|Acc]);
scaleaxis([],_,_,Acc) ->
   reverse(Acc).


pythag(_,0) ->
        0;
pythag(X,Y) ->
   case {abs(X),abs(Y)} of
           {A,B} when A>B -> A*math:sqrt(1+(B*B/A/A));
           {A,B} when A=<B -> B*math:sqrt(1+(A*A/B/B))
   end.
