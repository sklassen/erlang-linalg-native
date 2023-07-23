-module(linalg).
-vsn('1.2.0').
-author('simon.klassen').

-export([row/2, col/2, cell/3, set_cell/4, set_row/3, set_col/3]).
-export([transpose/1, t/1, flipud/1, fliplr/1]).
-export([det/1, inv/1, shape/1, reshape/2]).
-export([dot/2, inner/2, outer/2, matmul/2, solve/2]).
-export([zeros/1, ones/1, sequential/1, random/1]).
-export([zeros/2, ones/2, sequential/2, random/2,random/3]).
-export([fill/2, fill/3]).
-export([identity/1, diag/1, eye/1, eye/2]).
-export([add/2, sub/2, mul/2, divide/2, pow/2]).
-export([log/1,log10/1,log2/1, sqrt/1]).
-export([acos/1, asin/1, atan/1, acosh/1, asinh/1, atanh/1]).
-export([cos/1, cosh/1, sin/1, sinh/1, tan/1, tanh/1]).
-export([mean/1, median/1, std/1, var/1, cov/2,cor/2,cov2cor/1]).
-export([epsilon/1, exp/1, abs/1]).
-export([floor/1,ceil/1,around/1,around/2]).
-export([sum/1, sumsq/1, prod/1, norm/1]).
-export([roots/1, lu/1, qr/1, cholesky/1, svd/1]).
-export([real/1,imag/1]).
-export([polyfit/3, polyval/2]).
-export([min/1, max/1, argmin/1, argmax/1]).

-define(EPSILON, 1.0e-12).
-define(NA, na).
-define(ERR, err).

-type dim() :: non_neg_integer().
-type scalar() :: number().
-type vector() :: list(scalar()).
-type matrix() :: list(vector()).

% linalg shape
shape(X) when is_number(X) ->
    {};
shape([X | _] = Vector) when is_number(X) ->
    {length(Vector)};
shape([[X | _] | _] = Matrix) when is_number(X) ->
    NRows = length(Matrix),
    NCols = length(lists:nth(1, Matrix)),
    {NRows, NCols}.

-spec reshape(vector(),{dim(), dim()}) -> matrix().
reshape(Xs,{NR,NC})->
  reshape(lists:flatten(Xs),{NR,NC},[]).
reshape([],{_R,_C},Acc)->
  lists:reverse(Acc);
reshape(Xs,{NR,NC},Acc)->
  {Row,Rest}=lists:split(NC,Xs),
  reshape(Rest,{NR,NC},[Row|Acc]).

% generation (vector)
-spec zeros(dim()) -> vector().
zeros(0) ->
    [];
zeros(N) ->
    [0 || _ <- lists:seq(1, N)].

-spec ones(dim()) -> vector().
ones(0) ->
    [];
ones(N) ->
    [1 || _ <- lists:seq(1, N)].

-spec sequential(dim()) -> vector().
sequential(0) ->
    [];
sequential(N) ->
    [X || X <- lists:seq(1, N)].

-spec random(dim()) -> vector().
random(N) when is_integer(N) andalso N>0->
    random(N,rand:seed(exsplus),[]).
-spec random(dim(), dim()) -> matrix().
random(NR, NC) when is_integer(NR) andalso is_integer(NC) ->
    random(NR,NC,[{seed,rand:seed(exsplus)}]);

random(N,[{seed,Seed}]) when is_integer(N) ->
    random(N,Seed,[]).

random(NR,NC,[{seed,Seed}]) when is_integer(NR) andalso is_integer(NC) ->
    reshape(random(NR*NC,Seed,[]),{NR,NC});

random(0,_,Acc)->
    Acc;
random(N,Seed,Acc)->
    {X,NextSeed}=rand:uniform_s(Seed),
    random(N-1,NextSeed,[X|Acc]).


-spec fill(dim(), scalar()) -> vector().
fill(0, _) ->
    [];
fill(N, Value) ->
    [Value || _ <- lists:seq(1, N)].

% generation (matrix)
-spec zeros(dim(), dim()) -> matrix().
zeros(0, _) ->
    [[]];
zeros(_, 0) ->
    [[]];
zeros(NR, NC) ->
    [[0 || _ <- lists:seq(1, NC)] || _ <- lists:seq(1, NR)].

-spec ones(dim(), dim()) -> matrix().
ones(0, _) ->
    [[]];
ones(_, 0) ->
    [[]];
ones(NR, NC) ->
    [[1 || _ <- lists:seq(1, NC)] || _ <- lists:seq(1, NR)].

-spec sequential(dim(), dim()) -> matrix().
sequential(NR, NC) ->
    [[(((R - 1) * NC) + C) / 1.0 || C <- lists:seq(1, NC)] || R <- lists:seq(1, NR)].

-spec fill(dim(), dim(), scalar()) -> matrix().
fill(0, _, _) ->
    [[]];
fill(_, 0, _) ->
    [[]];
fill(NR, NC, Value) ->
    [fill(NC, Value) || _ <- lists:seq(1, NR)].

-spec eye(dim()) -> matrix().
eye(N) ->
    eye(N, N).

-spec eye(dim(), dim()) -> matrix().
eye(N, M) ->
    [
        [
            case {R, C} of
                {C, R} -> 1;
                {R, C} -> 0
            end
         || R <- lists:seq(1, M)
        ]
     || C <- lists:seq(1, N)
    ].

-spec diag(vector() | matrix()) -> matrix() | vector().
% (V)ector V->M
diag([X | _] = V) when is_number(X) ->
    [
        [
            case R==C of
                true -> lists:nth(R, V);
                false -> 0
            end
         || R <- lists:seq(1, length(V))
        ]
     || C <- lists:seq(1, length(V))
    ];
% (M)atrix M->V
diag([[X | _] | _] = M) when is_number(X) ->
    [lists:nth(R, lists:nth(R, M)) || R <- lists:seq(1, length(M))].

-spec identity(dim()) -> matrix().
identity(N) ->
    diag(ones(N)).

% Reference
-spec row(number(), matrix()) -> vector() | matrix().
row(0, _) ->
    [];
row(I, Matrix) when I > 0 ->
    lists:nth(I, Matrix);
row(I, Matrix) when I < 0 ->
    {A, [_ | B]} = lists:split(-(I + 1), Matrix),
    lists:append(A, B).

-spec col(number(), matrix()) -> vector() | matrix().
col(0, _) ->
    [];
col(J, Matrix) when J > 0 ->
    [lists:nth(J, Row) || Row <- Matrix];
col(J, Matrix) when J < 0 ->
    Colbind = fun({A, [_ | B]}) -> lists:append(A, B) end,
    [Colbind(lists:split(-(J + 1), Row)) || Row <- Matrix].

-spec cell(dim(), dim(), matrix()) -> scalar().
cell(I, J, Matrix) ->
    lists:nth(J, row(I, Matrix)).

-spec set_cell(dim(), dim(), scalar(), matrix()) -> matrix().
set_cell(I, J, Value, Matrix) ->
    Row = row(I, Matrix),
    set_nth(I, Matrix, set_nth(J, Row, Value)).

-spec set_row(dim(), vector() | matrix(), matrix()) -> matrix().
set_row(I, Value, Matrix) ->
    set_nth(I, Matrix, Value).

-spec set_col(dim(), vector() | matrix(), matrix()) -> matrix().
set_col(I, Value, Matrix) ->
    transpose(set_row(I, Value, transpose(Matrix))).

% Transformation 
-spec transpose(matrix()) -> matrix().
transpose([[]]) ->
    [];
transpose([[X]]) ->
    [[X]];
transpose([[] | Rows]) ->
    transpose(Rows);
transpose([[X | Xs] | Rows]) ->
    [[X | [H || [H | _] <- Rows]] | transpose([Xs | [Tail || [_ | Tail] <- Rows]])].

-spec t(matrix()) -> matrix().
t(M) ->
    transpose(M).

-spec flipud(matrix()) -> matrix().
flipud(M) ->
    lists:reverse(M).

-spec fliplr(matrix()) -> matrix().
fliplr(M) ->
    [lists:reverse(R) || R <- M].

% Sum Product (slower than inner, for big vectors, but succient)
-spec dot(vector(), vector()) -> scalar().
dot(VecA, VecB) ->
    lists:foldl(fun(X, Sum) -> Sum + X end, 0, lists:zipwith(fun(X, Y) -> X * Y end, VecA, VecB)).

% Matrix Multiplication
-spec matmul(matrix()|vector(), matrix()|vector()) -> matrix().
matmul(M = [V0|_], V = [X0|_]) when is_list(V0) andalso is_number(X0) ->
    [[ dot(Row,V) || Row <- M]];
matmul(V = [X0|_], M = [V0|_]) when is_number(X0) andalso is_list(V0) ->
    [[ dot(Row,V) || Row <- transpose(M)]];
matmul(M1 = [H1 | _], M2) when length(H1) =:= length(M2) ->
    matmul(M1, transpose(M2), []);
matmul([_H1 | _],_M2)->
  erlang:error("linalg:matmul matrix dim mismatch").
matmul([], _, R) ->
    lists:reverse(R);
matmul([Row | Rest], M2, R) ->
    matmul(Rest, M2, [outer(Row, M2) | R]).

-spec inner(vector(), vector()) -> scalar().
inner([], [], Sum) ->
    Sum;
inner([A], [B], Sum) ->
    Sum + A * B;
inner([A | VecA], [B | VecB], Sum) ->
    inner(VecA, VecB, Sum + A * B).
inner(VecA, VecB) ->
    inner(VecA, VecB, 0).

-spec outer(vector(), matrix()) -> vector().
outer(V1, V2) ->
    outer(V1, V2, []).
outer(_, [], R) ->
    lists:reverse(R);
outer(Row, [Col | Rest], R) ->
    outer(Row, Rest, [inner(Row, Col) | R]).

% Arithmetric
exp(M) ->
    sig1(M, fun(X) -> math:exp(X) end, []).

abs(M) ->
    sig1(M, fun(X) -> erlang:abs(X) end, []).

log(M) ->
    sig1(M, fun(X) -> math:log(X) end, []).

log10(M)->
    sig1(M, fun(X) -> math:log10(X) end, []).
log2(M)->
    sig1(M, fun(X) -> math:log2(X) end, []).

sqrt(M) ->
    sig1(M, fun(X) -> math:sqrt(X) end, []).

floor(M) ->
    sig1(M, fun(X) -> math:floor(X) end, []).

ceil(M) ->
    sig1(M, fun(X) -> math:ceil(X) end, []).

around(M) ->
    sig1(M, fun(X) -> round(X) end, []).

around(M,DP) ->
    Scale=math:pow(10,DP),
    sig1(M, fun(X) -> round(X*Scale)/Scale end, []).

epsilon(M) ->
    sig1(
        M,
        fun(X) ->
            case (erlang:abs(X) < ?EPSILON) of
                true -> 0;
                false -> X
            end
        end,
        []
    ).

% Trig 
acos(M)->
    sig1(M, fun(X) -> math:acos(X) end, []).
acosh(M)->
    sig1(M, fun(X) -> math:acosh(X) end, []).
asin(M)->
    sig1(M, fun(X) -> math:asin(X) end, []).
asinh(M)->
    sig1(M, fun(X) -> math:asinh(X) end, []).
atan(M)->
    sig1(M, fun(X) -> math:atan(X) end, []).
atanh(M)->
    sig1(M, fun(X) -> math:atanh(X) end, []).
cos(M)->
    sig1(M, fun(X) -> math:cos(X) end, []).
cosh(M)->
    sig1(M, fun(X) -> math:cosh(X) end, []).
sin(M)->
    sig1(M, fun(X) -> math:sin(X) end, []).
sinh(M)->
    sig1(M, fun(X) -> math:sinh(X) end, []).
tan(M)->
    sig1(M, fun(X) -> math:tan(X) end, []).
tanh(M)->
    sig1(M, fun(X) -> math:tanh(X) end, []).

real(M)->
    sig1(M, fun(X) -> case X of {R,_I}->R; _-> X end end, []).

imag(M)->
    sig1(M, fun(X) -> case X of {_R,I}->I; _-> 0 end end, []).

add(M1, M2) ->
    sig2(M1, M2, fun(A, B) -> A + B end, []).

sub(M1, M2) ->
    sig2(M1, M2, fun(A, B) -> A - B end, []).

mul(M1, M2) ->
    sig2(M1, M2, fun(A, B) -> A * B end, []).

divide(M1, M2) ->
    sig2(
        M1,
        M2,
        fun(A, B) ->
            case B==0 of
                true -> ?NA;
                false -> A / B
            end
        end,
        []
    ).

pow(M1, M2) ->
    sig2(M1, M2, fun(A, B) -> math:pow(A, B) end, []).


% Reductions
sum(X) when is_number(X) -> X;
sum(Xs) -> reduction(Xs, fun(X, Sum) -> Sum + X end, 0).

sumsq(X) when is_number(X) -> X * X;
sumsq(Xs) -> reduction(Xs, fun(X, SumSq) -> SumSq + X * X end, 0).

prod(X) when is_number(X) -> X;
prod(Xs) -> reduction(Xs, fun(X, Acc) -> Acc * X end, 1).

min(X) when is_number(X) -> X;
min([H | Vector]) when is_number(H) ->
    lists:foldl(fun(X, Min) -> erlang:min(X, Min) end, H, Vector);
min([H | Tail]) ->
    min(lists:flatten([H | Tail])).

max(X) when is_number(X) -> X;
max([H | Vector]) when is_number(H) ->
    lists:foldl(fun(X, Max) -> erlang:max(X, Max) end, H, Vector);
max([H | Tail]) ->
    max(lists:flatten([H | Tail])).

argmin(X) when is_number(X) -> 0;
argmin([H | Vector]) when is_number(H) ->
    IdxVector = lists:zip(lists:seq(2, length(Vector) + 1), Vector),
    element(
        1, lists:foldl(fun({I, X}, {K, Min}) -> argmin({I, X}, {K, Min}) end, {1, H}, IdxVector)
    );
argmin([H | Tail]) ->
    argmin(lists:flatten([H | Tail])).

argmin({I, X}, {_, Min}) when X < Min ->
    {I, X};
argmin({_, _}, {K, Min}) ->
    {K, Min}.

argmax(X) when is_number(X) -> 0;
argmax([H | Vector]) when is_number(H) ->
    IdxVector = lists:zip(lists:seq(2, length(Vector) + 1), Vector),
    element(
        1, lists:foldl(fun({I, X}, {K, Max}) -> argmax({I, X}, {K, Max}) end, {1, H}, IdxVector)
    );
argmax([H | Tail]) ->
    argmax(lists:flatten([H | Tail])).

argmax({I, X}, {_, Max}) when X > Max ->
    {I, X};
argmax({_, _}, {K, Max}) ->
    {K, Max}.

norm(X) when is_number(X) -> X;
norm(X) -> math:sqrt(sumsq(X)).

% Stats

mean(X) when is_number(X) -> X;
mean([X]) when is_number(X) -> X;
mean([X | Xs]) when is_number(X) ->
    sum([X | Xs]) / length([X | Xs]);
mean([[H | _] | _] = Matrix) when is_number(H) ->
    mean(lists:flatten(Matrix)).

median(X) when is_number(X) -> X;
median([X]) when is_number(X) -> X;
median([X | Xs]) when is_number(X) ->
    median_(lists:sort([X | Xs]));
median([[H | _] | _] = Matrix) when is_number(H) ->
    median(lists:flatten(Matrix)).

median_([X]) -> X;
median_([X1, X2]) -> (X1 + X2) / 2;
median_([_ | Tail]) -> median_(lists:droplast(Tail)).

std(X) ->
    math:sqrt(var(X)).

-spec var(vector()) -> scalar().
var(X) when is_number(X) -> 0;
var([X]) when is_number(X) -> 0;
var([X | Xs]) when is_number(X) ->
    {_Mean, SumSq, N} = welford(Xs, {X, 0, 1}),
    SumSq / N;
var([[H | _] | _] = Matrix) when is_number(H) ->
    var(lists:flatten(Matrix)).

% Welford
welford([], {Mean, SumSq, N}) ->
    {Mean, SumSq, N};
welford([H | Tail], {Mean, SumSq, N}) ->
    D1 = H - Mean,
    NewMean = Mean + D1 / (N + 1),
    D2 = H - NewMean,
    NewSumSq = SumSq + D1 * D2,
    welford(Tail, {NewMean, NewSumSq, N + 1}).

-spec cov(vector(), vector()) -> matrix().
cov([X], [Y]) when is_number(X) andalso is_number(Y) ->
    [[?NA, ?NA], [?NA, ?NA]];
cov([X | Xs], [Y | Ys]) when is_number(X) andalso is_number(Y) ->
    cov(Xs, Ys, {X, Y, 0, 0, 0, 0, 0, 1}).

cov([], [], {_X0, _Y0, DX, DXX, DY, DYY, DXY, N}) ->
    XX = (DXX - DX * DX / N) / (N - 1),
    YY = (DYY - DY * DY / N) / (N - 1),
    XY = (DXY - DX * DY / N) / (N - 1),
    [[XX, XY], [XY, YY]];
cov([X | XTail], [Y | YTail], {X0, Y0, DX, DXX, DY, DYY, DXY, N}) ->
    cov(
        XTail,
        YTail,
        {X0, Y0, 
         DX + (X - X0), DXX + (X - X0) * (X - X0), DY + (Y - Y0), 
         DYY + (Y - Y0) * (Y - Y0), DXY + (X - X0) * (Y - Y0), 
         N + 1}
    ).

-spec cor(vector(),vector())->matrix().
cor(Xs,Ys) ->
  pearson(Xs,Ys).

pearson(Xs,Ys) ->
  [[XX,XY],[XY,YY]]= cov(Xs,Ys),
  XY / (sqrt(XX) * sqrt(YY)).


-spec cov2cor(matrix()) -> matrix().
cov2cor(Cov) ->
  D = diag([1/X || X<- sqrt(diag(Cov))]),
  matmul(t(matmul(t(Cov),D)),D).

% Solves
-spec det(matrix()) -> scalar().
det(X) ->
    linalg_det:det(X).

-spec solve(matrix(), matrix()) -> matrix().
solve(X, B) ->
    Inv = inv(X),
    matmul(Inv, B).

-spec inv(matrix()) -> matrix().
inv([[X]]) ->
    [[1.0 / X]];
inv([[A, B], [C, D]]) ->
    case det([[A, B], [C, D]]) of
        0.0 -> ?ERR;
        Det -> [[D / Det, -B / Det], [-C / Det, A / Det]]
    end;
inv(Matrix) ->
    {NRows, NCols} = shape(Matrix),
    case det(Matrix) of
        0.0 -> ?ERR;
        Det -> divide(transpose(mul(minors(Matrix,NRows, NCols), cofactors(NRows, NCols))), Det)
    end.

-spec minors(matrix(),dim(),dim()) -> matrix().
minors(Matrix,NRows, NCols) ->
    [[det(col(-J, row(-I, Matrix))) || J <- lists:seq(1, NCols)] || I <- lists:seq(1, NRows)].

-spec cofactors(dim(),dim()) -> matrix().
cofactors(NRows, NCols) ->
    [[pow(-1, I) * pow(-1, J) || J <- lists:seq(0, NCols - 1)] || I <- lists:seq(0, NRows - 1)].

-spec roots(vector()) -> vector().
roots(Vector) ->
    linalg_roots:roots(Vector).

-spec cholesky(matrix()) -> matrix().
cholesky(RowWise) ->
    linalg_cholesky:cholesky(RowWise).

-spec lu(matrix()) -> {matrix(), matrix(), matrix()}.
lu(RowWise) ->
    linalg_lu:lu(RowWise).

-spec qr(matrix()) -> {matrix(), matrix()}.
qr(RowWise) ->
    linalg_qr:qr(RowWise).

-spec svd(matrix()) -> {matrix(), vector(),matrix()}.
svd(RowWise) ->
    linalg_svd_power:svd(RowWise).

-spec polyfit(vector(),vector(),integer()) -> vector().
polyfit(Xs,Ys,N) ->
    linalg_poly:polyfit(Xs,Ys,N).

-spec polyval(vector(),scalar()) -> scalar().
polyval(Coeff,X) ->
    linalg_poly:polyfit(Coeff,X).

% private arithmetic functions

sig1(X,_Fun, _) when is_atom(X) ->
    X;
sig1({R,I}, Fun, _) when is_number(R) andalso is_number(I) ->
    Fun({R,I});
sig1(X, Fun, _) when is_number(X) ->
    Fun(X);
sig1([], _Fun, Acc) ->
    lists:reverse(Acc);
sig1([H | _] = Vector, Fun, []) when is_number(H) orelse is_atom(H) orelse is_tuple(H)->
    [case is_atom(X) of true->X; false -> Fun(X) end || X <- Vector];
sig1([R1 | Matrix], Fun, Acc) ->
    sig1(Matrix, Fun, [[Fun(X) || X <- R1] | Acc]).

sig2([], [], _Fun, Acc) ->
    lists:reverse(Acc);
sig2(X, [], _Fun, Acc) when is_number(X) ->
    lists:reverse(Acc);
sig2([], X, _Fun, Acc) when is_number(X) ->
    lists:reverse(Acc);
sig2(A, B, Fun, []) when is_number(A) andalso is_number(B) ->
    Fun(A, B);
sig2(A, [B | Vector], Fun, []) when is_number(A) andalso is_number(B) ->
    [Fun(A, X) || X <- [B | Vector]];
sig2([A | Vector], B, Fun, []) when is_number(A) andalso is_number(B) ->
    [Fun(X, B) || X <- [A | Vector]];
sig2([A | VectorA], [B | VectorB], Fun, []) when is_number(A) andalso is_number(B) ->
    [Fun(X, Y) || {X, Y} <- lists:zip([A | VectorA], [B | VectorB])];
sig2(A, [R2 | M2], Fun, Acc) when is_number(A) ->
    sig2(A, M2, Fun, [[Fun(A, B) || B <- R2] | Acc]);
sig2([R1 | M1], B, Fun, Acc) when is_number(B) ->
    sig2(M1, B, Fun, [[Fun(A, B) || A <- R1] | Acc]);
sig2([R1 | M1], [R2 | M2], Fun, Acc) ->
    sig2(M1, M2, Fun, [[Fun(A, B) || {A, B} <- lists:zip(R1, R2)] | Acc]).

reduction([], _Fun, Init) ->
    Init;
reduction([X | _] = Vector, Fun, Init) when is_number(X) ->
    lists:foldl(Fun, Init, Vector);
reduction([V | _] = Matrix, Fun, Init) when is_list(V) ->
    reduction(lists:flatten(Matrix), Fun, Init).

set_nth(1, [_ | T], Value) -> [Value | T];
set_nth(I, [H | T], Value) -> [H | set_nth(I - 1, T, Value)].
