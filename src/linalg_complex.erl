-module(linalg_complex).
-vsn('1.0').
-author('tiago.klassen').
-define(SMALL, 1.0e-10).

-export([to_polar/1, from_polar/1]).
-export([exp/1, ln/1, pow/2, sqrt/1, qbrt/1, absolute/1, arg/1]). 
-export([add/1, add/2, mltp/1, mltp/2, reciprocal/1, usort/1]).

-type non_neg_real() :: number(). % but not negative please.
-type angle() :: number(). % principal value is -PI < argZ =< PI
-type complex_tuple() :: {Real::number(), Imaginary::number()}.
-type complex() :: Real::number() | complex_tuple().
-type polar() :: {polar, Radius::non_neg_real(), Angle::angle()}.

-spec from_polar(polar()) -> complex().
from_polar({polar, Radius, Angle}) ->
    to_complex({Radius*math:cos(Angle), Radius*math:sin(Angle)}).

-spec to_polar(complex()) -> polar().
to_polar(Complex) ->
    {polar, absolute(Complex), arg(Complex)}.

-spec exp(complex()) -> complex().
exp(Complex) ->
    {Real, Imaginary} = to_tuple(Complex),
    from_polar({polar, math:exp(Real), Imaginary}).

%% ln(X) aka log_e
-spec ln(complex()) -> complex().
ln(Complex) ->
    {polar, Radius, Angle} = to_polar(Complex),
    if
        Radius < ?SMALL -> throw(complex_infinity);
        true -> to_complex({math:log(Radius), Angle})
    end.

-spec pow(complex(), complex()) -> complex().
pow(Z, N) ->
    case {to_polar(Z), to_polar(N)} of
        {{_,Rz,_},{_,Rn,_}} when Rz<?SMALL,Rn<?SMALL ->
        	error({indeterminate_form, {pow, [Z,N]}});
        {{_,Rz,_},_} when Rz<?SMALL ->
            case to_tuple(N) of
                {Real,_} when Real>0 -> 0.0;
                _ -> throw(complex_infinity)
            end;
        _ -> exp(mltp(N, ln(Z)))
    end.

-spec sqrt(complex()) -> complex().
sqrt(Complex) ->
    pow(Complex, 1/2).

-spec qbrt(complex()) -> complex().
qbrt(Complex) ->
    pow(Complex, 1/3).


-spec absolute(complex()) -> number().
absolute({R,I}) when is_number(R), is_number(I) ->
    0.0+math:sqrt(R*R+I*I);
absolute(R) when is_number(R) ->
    0.0+abs(R).

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
    {polar, Ar, Aa} = to_polar(A),
    {polar, Br, Ba} = to_polar(B),
    from_polar({polar, Ar*Br, Aa+Ba}).

-spec reciprocal(complex()) -> complex().
reciprocal(A) ->
    {Real, Imaginary} = to_tuple(A),
    R2I2 = Real*Real + Imaginary*Imaginary,
    to_complex({Real/R2I2, -Imaginary/R2I2}).

%% function to find the angle.
-spec arg(complex()) -> number().
arg(Real) when is_number(Real) ->
    if
        Real>=0 -> 0.0;
        Real<0 -> math:pi()
    end;
arg({Real,Imaginary}) when abs(Imaginary) < ?SMALL->
    if
        Real>=0 -> 0.0;
        Real<0 -> math:pi()
    end;
arg({Real,Imaginary}) ->
    R = Real / math:sqrt(Real*Real + Imaginary*Imaginary),
    I = Imaginary / math:sqrt(Real*Real + Imaginary*Imaginary),
    if
        I>0 -> math:acos(R);
        I<0 -> -math:acos(R)
    end.

-spec usort([complex()]) -> complex().
usort([]) -> [];
usort(List) when is_list(List) ->
    [Head|Tail] = lists:usort(List),
    FoldFun = fun(X, {[Y|Prev], C}) ->
        case absolute(add(X,mltp(-1,Y))) < ?SMALL of
            true -> {[mltp(1/(C+1),add(X,mltp(C,Y)))|Prev], C+1};
            false -> {[X|[Y|Prev]], 1}
        end
    end,
    {Res, _C} = lists:foldr(FoldFun, {[Head],1}, lists:reverse(Tail)),
    lists:reverse(Res).

-spec to_tuple(complex()) -> complex_tuple().
to_tuple({R, I}) when is_number(R), is_number(I) -> {R, I};
to_tuple(Real) when is_number(Real) -> {Real, 0.0}.

-spec to_complex(complex()) -> complex().
to_complex({R, I}) when abs(R) < ?SMALL, abs(I) < ?SMALL -> 0.0;
to_complex({Real, Imaginary}) when abs(Imaginary) < ?SMALL ->
    to_float(Real);
to_complex({Real, Imaginary}) when abs(Real) < ?SMALL ->
    {0.0, to_float(Imaginary)};
to_complex({R,I}) -> {to_float(R), to_float(I)};
to_complex(N) -> to_float(N).

-spec to_float(number()) -> float().
to_float(Integer) when is_integer(Integer) -> 0.0+Integer;
to_float(Float) when is_float(Float) -> Float.
