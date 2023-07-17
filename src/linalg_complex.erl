-module(linalg_complex).
-vsn('1.0').
-author('tiago.klassen').
-define(SMALL, 1.0e-10).

-export([sqrt/1, qbrt/1]). 
-export([add/1, add/2, mltp/1, mltp/2, reciprocal/1, min/1, max/1]).

-type complex_tuple() :: {Real::number(), Imaginary::number()}.
-type complex() :: Real::number() | complex_tuple().

-spec sqrt(complex()) -> complex().
sqrt(Real) when is_number(Real) ->
    sqrt({Real, 0});
sqrt({Real, Imaginary}) when abs(Imaginary) < ?SMALL ->
    if
        abs(Real) < ?SMALL -> 0;
        Real > 0 -> math:sqrt(Real);
        Real < 0 -> {0, math:sqrt(-Real)}
    end;
sqrt(Complex) ->
    % yes, there is a way to solve this. but not now.
    error(non_real_number_sqrt, Complex).

-spec qbrt(complex()) -> complex().
qbrt(Complex) ->
    W = {-1/2, sqrt(3)/2},
    Ans0 = qbrt1(Complex),
    Ans1 = mltp(Ans0, W),
    Ans2 = mltp(Ans1, W),
    if
        is_number(Ans0) -> Ans0;
        is_number(Ans1) -> Ans1;
        is_number(Ans2) -> Ans2;
        true -> max([Ans0,Ans1,Ans2])
    end.
qbrt1(Real) when is_number(Real) ->
    qbrt1({Real, 0});
qbrt1({Real, Imaginary}) when abs(Imaginary) < ?SMALL ->
    if
        abs(Real) < ?SMALL ->
            0;
        Real > 0 -> 
            math:pow(Real, 1/3);
        Real < 0 ->
            P = math:pow(-Real, 1/3),
            Q = math:pi()/3,
            to_complex({P*math:cos(Q), P*math:sin(Q)})
    end;
qbrt1({R, I}) when abs(R) < ?SMALL ->
    P = math:pow(abs(I), 1/3),
    Q = I/abs(I) * math:pi()/6,
    to_complex(mltp(P, {math:cos(Q), math:sin(Q)}));
qbrt1({R, I}) ->
    P = math:pow(R*R+I*I, 1/6),
    {Np, Q} = case {R>0, I>0} of
        {true, true} -> {1, 1/3*math:atan(I/R)};
        {false, true} -> {1, 1/3*(math:pi()-math:atan(I/abs(R)))};
        {true, false} -> {-1, 1/3*math:atan(abs(I)/R)};
        {false, false} -> {1, 1/3*(math:atan(abs(I)/abs(R))-math:pi())}
    end,
    to_complex({P*math:cos(Q), Np*P*math:sin(Q)}).

-spec add([complex(), ...]) -> complex().
add(ComplexList) ->
    FoldFun = fun(X, Prev) -> add(X, Prev) end,
    lists:foldl(FoldFun, 0, ComplexList).

-spec add(complex(), complex()) -> complex().
add(A, B) ->
    {Ar, Ai} = to_tuple(A),
    {Br, Bi} = to_tuple(B),
    to_complex({Ar+Br, Ai+Bi}).

-spec mltp([complex(), ...]) -> complex().
mltp(ComplexList) ->
    FoldFun = fun(X, Prev) -> mltp(X, Prev) end,
    lists:foldl(FoldFun, 1, ComplexList).

-spec mltp(complex(), complex()) -> complex().
mltp(A, B) ->
    {Ar, Ai} = to_tuple(A),
    {Br, Bi} = to_tuple(B),
    to_complex({Ar*Br-Ai*Bi, Ar*Bi+Br*Ai}).

-spec reciprocal(complex()) -> complex().
reciprocal(A) ->
    {Real, Imaginary} = to_tuple(A),
    R2I2 = Real*Real + Imaginary*Imaginary,
    to_complex({Real/R2I2, -Imaginary/R2I2}).

-spec angle(complex()) -> number().
angle(Real) when is_number(Real) ->
    if
        Real>0 -> 0;
        Real<0 -> math:pi()
    end;
angle({Real,Imaginary}) ->
    R = Real / sqrt(Real*Real + Imaginary*Imaginary),
    I = Imaginary / sqrt(Real*Real + Imaginary*Imaginary),
    if
        I>=0 -> math:acos(R);
        I<0 -> 2*math:pi() - math:acos(R)
    end.

-spec min([complex(), ...]) -> complex().
min([H|T]) ->
    to_complex(find_minmax(min, T, H)).

-spec max([complex(), ...]) -> complex().
max([H|T]) ->
    to_complex(find_minmax(max, T, H)).

find_minmax(MinMax, [], Ans) ->
    Ans;
find_minmax(MinMax, [Complex|Next], Prev) ->
    {C, Ci} = to_tuple(Complex),
    {P, Pi} = to_tuple(Prev),
    if
        MinMax==max, C > P -> find_minmax(MinMax, Next, Complex);
        MinMax==min, C < P -> find_minmax(MinMax, Next, Complex);
        MinMax==max, C==P, Ci > Pi -> find_minmax(MinMax, Next, Complex);
        MinMax==min, C==P, Ci < Pi -> find_minmax(MinMax, Next, Complex);
        true -> find_minmax(MinMax, Next, Prev)
    end.

-spec to_tuple(complex()) -> complex_tuple().
to_tuple({Real, Imaginary}) -> {Real, Imaginary};
to_tuple(Real) -> {Real, 0}.

-spec to_complex(complex()) -> complex().
to_complex({Real, Imaginary}) when abs(Imaginary) < ?SMALL -> Real;
to_complex(Complex) -> Complex.
