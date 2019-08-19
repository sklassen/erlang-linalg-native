-module(linalg_native). 
-vsn('1.0').
-author('simon.klassen').
-export([dot/2]). 
-export([transpose/1,matmul/2]). 

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
	SumProd = dot(lists:nth(R,M1),lists:nth(C,M2)),
	rowmult(I, C-1, R, [SumProd|L], M1, M2).

matmul(_, _, 0, MM, _, _) -> MM;
matmul(I, C, R, MM, M1, M2) ->
	NewRow = rowmult(I, C, R, [], M1, M2),
	matmul(I, C, R-1, [NewRow|MM], M1, M2).

matmul(M1, M2) -> 
	Inner = length(M2),
	NCols = length(lists:nth(1,M2)), 
	NRows = length(M1), 
	matmul(Inner, NCols, NRows,[], M1, transpose(M2)).

