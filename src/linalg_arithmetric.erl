-module(linalg_arithmetric). 
-vsn('1.0').
-author('simon.klassen').
-import(lists,[reverse/1,zip/2,foldl/3]).
-export([sum/1,norm/1]).
-export([epsilon/1,exp/1,log/1]).
-export([add/2,sub/2,mul/2,divide/2,pow/2]).
-define(EPSILON,1/1000000).

% Arithmetric
-spec exp(list())->list().
exp(M)->
    sig1(M,fun(X)->math:exp(X) end,[]).

log(M)->
    sig1(M,fun(X)->math:exp(X) end,[]).

epsilon(M)->
    sig1(M,fun(X)-> case (abs(X)<?EPSILON) of true->0; false->X end end,[]).

add(M1,M2)->
    sig2(M1,M2,fun(A,B)->A+B end,[]).

sub(M1,M2)->
    sig2(M1,M2,fun(A,B)->A-B end,[]).

mul(M1,M2)->
    sig2(M1,M2,fun(A,B)->A*B end,[]).

divide(M1,M2)->
    sig2(M1,M2,fun(A,B)->case B of 0-> na; B-> A/B end end,[]).

pow(M1,M2)->
    sig2(M1,M2,fun(A,B)->math:pow(A,B) end,[]).

% Reductions

norm([H|_]=Vector) when is_number(H)->
    math:sqrt(sum(pow(Vector,2))).

sum(X) when is_number(X)->X;
sum([])->0;
sum([H|_]=Vector) when is_number(H)->
    foldl(fun(A,Sum)->Sum+A end,0,Vector);
sum([H|Tail])->
    sum(H)+sum(Tail).

% private functions

sig1(X,Fun,_) when is_number(X)->
    Fun(X);
sig1([],_Fun,Acc)->
    reverse(Acc);
sig1([H|_]=Vector,Fun,[]) when is_number(H)->
    [Fun(X)||X <-Vector];
sig1([R1|Matrix],Fun,Acc)->
    sig1(Matrix,Fun,[[Fun(X)||X <-R1]|Acc]).


sig2([],[],_Fun,Acc)->
    reverse(Acc);
sig2(X,[],_Fun,Acc) when is_number(X)->
    reverse(Acc);
sig2([],X,_Fun,Acc) when is_number(X)->
    reverse(Acc);
sig2(A,B,Fun,[]) when is_number(A) andalso is_number(B)->
    Fun(A,B);
sig2(A,[B|Vector],Fun,[]) when is_number(A) andalso is_number(B)->
    [Fun(A,X)||X<-[B|Vector]];
sig2([A|Vector],B,Fun,[]) when is_number(A) andalso is_number(B)->
    [Fun(X,B)||X<-[A|Vector]];
sig2([A|VectorA],[B|VectorB],Fun,[]) when is_number(A) andalso is_number(B)->
    [Fun(X,Y)||{X,Y}<-zip([A|VectorA],[B|VectorB])];
sig2(A,[R2|M2],Fun,Acc) when is_number(A)->
    sig2(A,M2,Fun,[[Fun(A,B)||B<-R2]|Acc]);
sig2([R1|M1],B,Fun,Acc) when is_number(B)->
    sig2(M1,B,Fun,[[Fun(A,B)||A<-R1]|Acc]);
sig2([R1|M1],[R2|M2],Fun,Acc)->
    sig2(M1,M2,Fun,[[Fun(A,B)||{A,B}<-zip(R1,R2)]|Acc]).
