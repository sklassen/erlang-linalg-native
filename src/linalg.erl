-module(linalg).
-vsn('1.0.1').
-author('simon.klassen').

-import(lists, [reverse/1, append/2, nth/2, seq/2, split/2, zip/2, foldl/3]).

-export([row/2, col/2, cell/3, set_cell/4]).
-export([transpose/1, flipud/1, fliplr/1]).
-export([det/1, inv/1, shape/1]).
-export([dot/2, inner/2, outer/2, matmul/2, solve/2]).
-export([zeros/1, ones/1, sequential/1, random/1]).
-export([zeros/2, ones/2, sequential/2, random/2]).
-export([fill/2, fill/3]).
-export([identity/1, diag/1, eye/1, eye/2]).
-export([add/2, sub/2, mul/2, divide/2, pow/2]).
-export([mean/1, median/1, std/1, var/1, cov/2]).
-export([epsilon/1, exp/1, abs/1, log/1, sqrt/1]).
-export([sum/1, sumsq/1, prod/1, norm/1]).
-export([roots/1, qr/1, cholesky/1]).
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
    NCols = length(nth(1, Matrix)),
    {NRows, NCols}.

% generation (vector)
-spec zeros(dim()) -> vector().
zeros(0) ->
    [];
zeros(N) ->
    [0 || _ <- seq(1, N)].

-spec ones(dim()) -> vector().
ones(0) ->
    [];
ones(N) ->
    [1 || _ <- seq(1, N)].

-spec sequential(dim()) -> vector().
sequential(0) ->
    [];
sequential(N) ->
    [X || X <- seq(1, N)].

-spec random(dim()) -> vector().
random(0) ->
    [];
random(N) ->
    [rand:uniform() || _ <- seq(1, N)].

-spec fill(dim(), scalar()) -> vector().
fill(0, _) ->
    [];
fill(N, Value) ->
    [Value || _ <- seq(1, N)].

% generation (matrix)
-spec zeros(dim(), dim()) -> matrix().
zeros(0, _) ->
    [[]];
zeros(_, 0) ->
    [[]];
zeros(NR, NC) ->
    [[0 || _ <- seq(1, NC)] || _ <- seq(1, NR)].

-spec ones(dim(), dim()) -> matrix().
ones(0, _) ->
    [[]];
ones(_, 0) ->
    [[]];
ones(NR, NC) ->
    [[1 || _ <- seq(1, NC)] || _ <- seq(1, NR)].

-spec sequential(dim(), dim()) -> matrix().
sequential(NR, NC) ->
    [[(((R - 1) * NC) + C) / 1.0 || C <- seq(1, NC)] || R <- seq(1, NR)].

-spec random(dim(), dim()) -> matrix().
random(NR, NC) ->
    [[rand:uniform() || _ <- seq(1, NC)] || _ <- seq(1, NR)].

-spec fill(dim(), dim(), scalar()) -> matrix().
fill(0, _, _) ->
    [[]];
fill(_, 0, _) ->
    [[]];
fill(NR, NC, Value) ->
    [fill(NC, Value) || _ <- seq(1, NR)].

-spec eye(dim()) -> matrix().
eye(0) ->
    eye([[]]);
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
         || R <- seq(1, M)
        ]
     || C <- seq(1, N)
    ].

-spec diag(vector() | matrix()) -> matrix() | vector().
% (V)ector V->M
diag([X | _] = V) when is_number(X) ->
    [
        [
            case R of
                C -> nth(R, V);
                _ -> 0
            end
         || R <- seq(1, length(V))
        ]
     || C <- seq(1, length(V))
    ];
% (M)atrix M->V
diag([[X | _] | _] = M) when is_number(X) ->
    [nth(R, nth(R, M)) || R <- seq(1, length(M))].

-spec identity(dim()) -> matrix().
identity(N) ->
    diag(ones(N)).

% Reference
-spec row(dim(), matrix()) -> vector() | matrix().
row(0, _) ->
    [];
row(I, Matrix) when I > 0 ->
    nth(I, Matrix);
row(I, Matrix) when I < 0 ->
    {A, [_ | B]} = split(-(I + 1), Matrix),
    append(A, B).

-spec col(dim(), matrix()) -> vector() | matrix().
col(0, _) ->
    [];
col(J, Matrix) when J > 0 ->
    [nth(J, Row) || Row <- Matrix];
col(J, Matrix) when J < 0 ->
    Colbind = fun({A, [_ | B]}) -> append(A, B) end,
    [Colbind(split(-(J + 1), Row)) || Row <- Matrix].

-spec cell(dim(), dim(), matrix()) -> scalar().
cell(I, J, Matrix) ->
    nth(J, row(I, Matrix)).

-spec set_cell(dim(), dim(), scalar(), matrix()) -> matrix().
set_cell(I, J, Value, Matrix) ->
    Row = row(I, Matrix),
    set_nth(I, Matrix, set_nth(J, Row, Value)).

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

-spec flipud(matrix()) -> matrix().
flipud(M) ->
    reverse(M).

-spec fliplr(matrix()) -> matrix().
fliplr(M) ->
    [reverse(R) || R <- M].

% Sum Product (slower than inner, for big vectors, but succient)
-spec dot(vector(), vector()) -> scalar().
dot(VecA, VecB) ->
    foldl(fun(X, Sum) -> Sum + X end, 0, lists:zipwith(fun(X, Y) -> X * Y end, VecA, VecB)).

% Matrix Multiplication
-spec matmul(matrix(), matrix()) -> matrix().
matmul(M1 = [H1 | _], M2) when length(H1) =:= length(M2) ->
    matmul(M1, transpose(M2), []).
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

-spec outer(vector(), vector()) -> matrix().
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

sqrt(M) ->
    sig1(M, fun(X) -> math:sqrt(X) end, []).

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
            case B of
                0 -> ?NA;
                B -> A / B
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
        {X0, Y0, DX + (X - X0), DXX + (X - X0) * (X - X0), DY + (Y - Y0), DYY + (Y - Y0) * (Y - Y0),
            DXY + (X - X0) * (Y - Y0), N + 1}
    ).

% Solves
-spec det(matrix()) -> scalar().
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
% laplace
det([H | Tail]) ->
    sum([pow(-1, J - 1) * X * det(col(-J, Tail)) || {J, X} <- zip(seq(1, length(H)), H)]).

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
        Det -> [[D / Det, -1 / Det * B], [-1 / Det * C, A / Det]]
    end;
inv(M) ->
    case det(M) of
        0.0 -> ?ERR;
        Det -> divide(transpose(mul(minors(M), cofactors(M))), Det)
    end.

minors(Matrix) ->
    {NRows, NCols} = shape(Matrix),
    [[det(col(-J, row(-I, Matrix))) || J <- seq(1, NCols)] || I <- seq(1, NRows)].

cofactors(Matrix) ->
    {NRows, NCols} = shape(Matrix),
    [[pow(-1, I) * pow(-1, J) || J <- seq(0, NCols - 1)] || I <- seq(0, NRows - 1)].

-spec roots(vector()) -> vector().
roots(Vector) ->
    linalg_roots:roots(Vector).

-spec cholesky(matrix()) -> matrix().
cholesky(RowWise) ->
    linalg_cholesky:crout(RowWise).

-spec qr(matrix()) -> {matrix(), matrix()}.
qr(RowWise) ->
    linalg_qr:qr(RowWise).

% private arithmetic functions

sig1(X, Fun, _) when is_number(X) ->
    Fun(X);
sig1([], _Fun, Acc) ->
    reverse(Acc);
sig1([H | _] = Vector, Fun, []) when is_number(H) ->
    [Fun(X) || X <- Vector];
sig1([R1 | Matrix], Fun, Acc) ->
    sig1(Matrix, Fun, [[Fun(X) || X <- R1] | Acc]).

sig2([], [], _Fun, Acc) ->
    reverse(Acc);
sig2(X, [], _Fun, Acc) when is_number(X) ->
    reverse(Acc);
sig2([], X, _Fun, Acc) when is_number(X) ->
    reverse(Acc);
sig2(A, B, Fun, []) when is_number(A) andalso is_number(B) ->
    Fun(A, B);
sig2(A, [B | Vector], Fun, []) when is_number(A) andalso is_number(B) ->
    [Fun(A, X) || X <- [B | Vector]];
sig2([A | Vector], B, Fun, []) when is_number(A) andalso is_number(B) ->
    [Fun(X, B) || X <- [A | Vector]];
sig2([A | VectorA], [B | VectorB], Fun, []) when is_number(A) andalso is_number(B) ->
    [Fun(X, Y) || {X, Y} <- zip([A | VectorA], [B | VectorB])];
sig2(A, [R2 | M2], Fun, Acc) when is_number(A) ->
    sig2(A, M2, Fun, [[Fun(A, B) || B <- R2] | Acc]);
sig2([R1 | M1], B, Fun, Acc) when is_number(B) ->
    sig2(M1, B, Fun, [[Fun(A, B) || A <- R1] | Acc]);
sig2([R1 | M1], [R2 | M2], Fun, Acc) ->
    sig2(M1, M2, Fun, [[Fun(A, B) || {A, B} <- zip(R1, R2)] | Acc]).

reduction([], _Fun, Init) ->
    Init;
reduction([X | _] = Vector, Fun, Init) when is_number(X) ->
    foldl(Fun, Init, Vector);
reduction([V | _] = Matrix, Fun, Init) when is_list(V) ->
    reduction(lists:flatten(Matrix), Fun, Init).

set_nth(1, [_ | T], Value) -> [Value | T];
set_nth(I, [H | T], Value) -> [H | set_nth(I - 1, T, Value)].
