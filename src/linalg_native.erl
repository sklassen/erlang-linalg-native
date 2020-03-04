-module(linalg_native). 
-vsn('1.0').
-author('simon.klassen').
-import(math,[pow/2]).
-import(lists,[reverse/1,append/2,nth/2,seq/2,split/2,zip/2,foldl/3]).
-export([row/2,col/2]). 
-export([transpose/1,det/1,inv/1,shape/1,dot/2,matmul/2]). 

shape(Matrix)->
   NCols = length(nth(1,Matrix)), 
   NRows = length(Matrix), 
   {NRows,NCols}.

transpose([[]]) -> [];
transpose([[X]]) -> [[X]];
transpose([[] | XXs]) -> transpose(XXs);
transpose([[X | Xs] | XXs]) -> [[X | [H || [H | _Tail ] <- XXs]] | transpose([Xs | [Tail || [_|Tail] <- XXs]])].


dot([],[],_) ->[];
dot([A],[B],Sum) ->Sum+A*B;
dot([A|VecA],[B|VecB],Sum) ->dot(VecA,VecB,Sum+A*B).
dot(VecA,VecB) ->dot(VecA,VecB,0).

rowmult(_, 0, _, L, _, _) -> L;
rowmult(I, C, R, L, M1, M2) -> 
	SumProd = dot(nth(R,M1),nth(C,M2)),
	rowmult(I, C-1, R, [SumProd|L], M1, M2).

matmul(_, _, 0, MM, _, _) -> MM;
matmul(I, C, R, MM, M1, M2) ->
	NewRow = rowmult(I, C, R, [], M1, M2),
	matmul(I, C, R-1, [NewRow|MM], M1, M2).

matmul(M1, M2) -> 
	Inner = length(M2),
	NCols = length(nth(1,M2)), 
	NRows = length(M1), 
	matmul(Inner, NCols, NRows,[], M1, transpose(M2)).

row(I,Matrix) when I>0 ->
    [nth(I,Matrix)];
row(I,Matrix) when I<0 ->
    {A,[_|B]}=split(-(I+1),Matrix),
    append(A,B).

col(J,Matrix) ->
    transpose(row(J,transpose(Matrix))).

det([[X]])->
    X;
det([[A,B],[C,D]])->
    A*D-B*C;
det([H|Tail])->
    foldl(fun(A,Sum)->Sum+A end,0,[pow(-1,J-1)*X*det(col(-J,Tail))||{J,X}<-zip(seq(1,length(H)),H)]).

inv([[X]])->
    [[1.0/X]];
inv([[A,B],[C,D]])->
    case det([[A,B],[C,D]]) of
       0.0->err;
       Det->[[D/Det,-1/Det*B],[-1/Det*C,A/Det]]
    end;
inv(M)->
    case det(M) of
       0.0->err;
       Det->mul(transpose(mul(minors(M),cofactors(M))),[[1/Det]])
    end.

minors(Matrix)->
    {NRows,NCols}=shape(Matrix),
    [[det(col(-J,row(-I,Matrix)))||J<-seq(1,NCols)]||I<-seq(1,NRows)].

cofactors(Matrix)->
    {NRows,NCols}=shape(Matrix),
    [[pow(-1,I)*pow(-1,J)||J<-seq(0,NCols-1)]||I<-seq(0,NRows-1)].


% Arithmetric

mul(M1,M2)->
    mul(M1,M2,[]).
mul([],[],Acc)->
    reverse(Acc);
mul([[_]],[],Acc)->
    reverse(Acc);
mul([],[[_]],Acc)->
    reverse(Acc);
mul([[A]],[R2|M2],Acc)->
    mul([[A]],M2,[[A*B||B<-R2]|Acc]);
mul([R1|M1],[[B]],Acc)->
    mul(M1,[[B]],[[A*B||A<-R1]|Acc]);
mul([R1|M1],[R2|M2],Acc)->
    mul(M1,M2,[[A*B||{A,B}<-zip(R1,R2)]|Acc]).
