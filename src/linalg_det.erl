-module(linalg_det).
-export([det/1]).

% 1x1 
det([[X]]) ->
    X;
% well known 2x2 ...
det([[A, B], [C, D]]) ->
    A * D - B * C;
% ... rule of sarrus 3x3 ...
det([[A, B, C], [D, E, F], [G, H, I]]) ->
    A * E * I + B * F * G + C * D * H - C * E * G - B * D * I - A * F * H;
% ... and extention 4x4 (speeds up processing for larger matrix by 2-3x)
det([[A, B, C, D], [E, F, G, H], [I, J, K, L], [M, N, O, P]]) ->
    A * F * K * P - A * F * L * O - A * G * J * P + A * G * L * N + A * H * J * O - A * H * K * N -
    B * E * K * P + B * E * L * O +
    B * G * I * P - B * G * L * M - B * H * I * O + B * H * K * M + C * E * J * P -
    C * E * L * N - C * F * I * P + C * F * L * M +
    C * H * I * N - C * H * J * M - D * E * J * O + D * E * K * N + D * F * I * O -
    D * F * K * M - D * G * I * N + D * G * J * M;
% Laplace 5x5 
det([H | Tail]) when length(H)==5 andalso length([H|Tail])==5 ->
		linalg:sum([math:pow(-1, J - 1) * X * det(linalg:col(-J, Tail)) || {J, X} <- lists:zip(lists:seq(1, length(H)), H)]);
% Remaining square matrix
det([H | Tail]=RowWise) when length(H)==length([H|Tail]) ->
	det(RowWise,1);

det(_) ->
		erlang:error({error,matrix_not_square}).

% Gauss-Jordan elimination
det(RowWise,N) when length(RowWise)==N->
	linalg:prod(linalg:diag(RowWise));
det(RowWise,N)->
  {Top,Bottom} = lists:split(N, RowWise),
	Row = lists:last(Top),
	case lists:nth(N,Row) of
	  0 -> case swap(RowWise,N,N+1) of
					na -> 0;
					New-> -det(New,N)
				 end;
		H -> Col=linalg:divide(linalg:col(N,Bottom),H),
				 Sub = [ linalg:mul(Row,X) || X<-Col ],
				 det(lists:append(Top,linalg:sub(Bottom,Sub)),N+1)
	end.

swap(RowWise, _, J) when J > length(RowWise) ->
  na;
swap(RowWise, I, J) ->
	case linalg:cell(J,I,RowWise) of
		0 -> swap(RowWise, I, J+1);
		_ -> {A, [C | Tail]} = lists:split(I - 1, RowWise),
    		 {B, D} = lists:split(J - I, Tail),
    		 lists:append([A, B, [C], D])
	end.
