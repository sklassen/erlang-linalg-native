-module(linalg_minimize).
-define(GLD,0.618033989).
-define(SMALL,0.00001).
-export([minimize/2,minimize/3]).

minimize(Fun,Para)->
  minimize(Fun,Para,golden).

minimize(Fun,Para,Options)->
  case Options of
    golden->golden(Fun,Para);
    brute->brute(Fun,Para)
  end.

golden(Fun,{A,B}) when A>B ->
  golden(Fun,{B,A});
golden(_Fun,{A,B}) when abs(A-B)<?SMALL ->
  (A+B)/2;
golden(Fun,{A,B})->
  X1=?GLD*A+(1.0-?GLD)*B,
  X2=?GLD*B+(1.0-?GLD)*A,
  try {Fun(X1),Fun(X2)} of
    {FX1,FX2} when FX1<FX2 -> golden(Fun,{A,X2});
    {_,_} -> golden(Fun,{X1,B})
  catch 
    _:_ -> erlang:error("minimize function failed (golden)")
  end.

brute(Fun,{A,B}) when A>B ->
  brute(Fun,{B,A});
brute(Fun,{A,B}) ->
  N=1000,
  Step=max(1,round((B/?SMALL-A/?SMALL)/N)),
  Steps = lists:seq(round(min(A,B)/?SMALL),round(max(A,B)/?SMALL),Step),
  brute(Fun,Steps,{A,Fun(A)}).
brute(_Fun,[],{ArgMin,_Min}) ->
  ArgMin;
brute(Fun,[X|Tail],{ArgMin,Min})->
  try Fun(X*?SMALL) of
    NewMin when NewMin<Min -> brute(Fun,Tail,{X*?SMALL,NewMin});
    _Next -> brute(Fun,Tail,{ArgMin,Min})
  catch 
    _:_ -> erlang:error("minimize function failed (golden)")
  end.
