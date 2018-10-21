%% Doesn't check matrix lengths match but 
%% not restricted to square matrices.
%%
%% Represent matrix as tuple of tuples.
%% Doesn't use destructive assignment so building a
%% matrix creates an awful lot of temporary lists.
%%
%% Usage: start from command line with
%%     erlc mmultiply.erl
%%     erl -noinput -s mmultiply main 1

-module(matrix). 
-vsn('2011.1201').
-author('simon.klassen@adyne.com').
-export([zero/2,sequential/2,identity/1,diag/1,transpose/1,mmultiply/2]). 
-import(lists, [reverse/1,nth/2,seq/2]).

transpose([]) -> [];
transpose([[] | Xss]) -> transpose(Xss);
transpose([[X | Xs] | Xss]) -> [[X | [H || [H | _T ] <- Xss]] | transpose([Xs | [T || [_H | T] <- Xss]])].

sumprod([],[],_) ->[];
sumprod([A],[B],Sum) ->Sum+A*B;
sumprod([A|VecA],[B|VecB],Sum) ->sumprod(VecA,VecB,Sum+A*B).
sumprod(VecA,VecB) ->sumprod(VecA,VecB,0).

rowmult(_, 0, _, L, _, _) -> L;
rowmult(I, C, R, L, M1, M2) -> 
    SumProd = sumprod(nth(R,M1),nth(C,M2)),
    rowmult(I, C-1, R, [SumProd|L], M1, M2).

mmultiply(_, _, 0, MM, _, _) -> MM;
mmultiply(I, C, R, MM, M1, M2) ->
    NewRow = rowmult(I, C, R, [], M1, M2),
    mmultiply(I, C, R-1, [NewRow|MM], M1, M2).

mmultiply(M1, M2) -> 
    Inner = length(M2),
    NCols = length(nth(1,M2)), 
    NRows = length(M1), 
    mmultiply(Inner, NCols, NRows,[], M1, transpose(M2)).

zero(NR,NC) -> 
	[ [ 0.0 ||_<-seq(1,NC)] || _<-seq(1,NR)].

sequential(NR,NC) -> 
	[ [ (((R-1)*NC)+C)/1.0 || C<-seq(1,NC)] || R<-seq(1,NR)].

identity(N) -> 
	[ [ case R of C -> 1; _->0 end||R<-seq(1,N)] || C<-seq(1,N)].

diag(X) -> 
	[ [ case R of C -> nth(R,X); _->0.0 end||R<-seq(1,length(X))] || C<-seq(1,length(X))].



