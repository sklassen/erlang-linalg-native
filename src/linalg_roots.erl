-module(linalg_roots).
-vsn('1.0').
-author('simon.klassen').
-import(linalg_complex, [sqrt/1, qbrt/1, absolute/1, usort/1]).
-import(linalg_complex, [add/1, add/2, mltp/1, mltp/2, reciprocal/1]).
-export([roots/1]).
-define(SMALL, 1.0e-10).

-type complex() :: Real::number() | {Real::number(), Imaginary::number()}.
-type vector() :: list(complex()).
-spec roots(vector()) -> vector().
% null case
roots([]) ->
    [];
% singleton
roots([_]) ->
    [];
% ignore higher order when alpha is small
roots([{R,I} | Tail]) when abs(R) < ?SMALL, abs(I) < ?SMALL ->
    roots(Tail);
roots([A | Tail]) when is_number(A), abs(A) < ?SMALL ->
    roots(Tail);
% ax + b = 0
roots([A, B]) ->
    [mltp([-1, B, reciprocal(A)])];
% quadratic polynomial
% well known closed form
% ax^2 + bx + c = 0
roots([A, B, C]) ->
    Denom = reciprocal(mltp(2,A)), % 1/(2*A)
    P = mltp([-1, B, Denom]), % -b/(2*A)
    % Q = sqrt(B^2-4*A*C)/(2*A)
    Q = mltp(Denom, sqrt(add(mltp(B, B), mltp([-4, A, C])))),
    usort([add(P, Q), add(P, mltp(-1, Q))]); % P +/- Q
% cubic polynomial
% viete
% x^3 + bx^2 + cx +d = 0
roots([1, B, C, D]) ->
    P = mltp(-1/3, B), % -b/3
    % Q = (-2B^3+9BC-27D)/54
    Q = mltp(1/54, add([mltp([-2,B,B,B]), mltp([9, B, C]), mltp(-27, D)])),
    % R = sqrt(3*(27D^2-18BCD+4B^3D+4C^3-B^2C^2))/18
    R = mltp(1/18, sqrt(mltp(3, add([
        mltp([27,D,D]), mltp([-18,B,C,D]),
        mltp([4,B,B,B,D]), mltp([4,C,C,C]), mltp([-1,B,B,C,C])
    ])))),
    
    % W = (-1+sqrt(3)i)/2
    W = {-1/2, sqrt(3)/2},
    
    % QpR = qbrt(Q+R)
    % QnR = qbrt(Q-R)
    QpR0 = qbrt(add(Q, R)),
    QnR0 = qbrt(add(Q, mltp(-1, R))),
    QpRs = [QpR0, mltp(W,QpR0), mltp([W,W,QpR0])],
    QnRs = [QnR0, mltp(W,QnR0), mltp([W,W,QnR0])],
    
    % X = P + QpRs + QnRs
    PossibleX = [add([P, QpR, QnR]) || QpR<-QpRs, QnR<-QnRs],
    XList = lists:filter(roots_filter([1,B,C,D]), PossibleX),
    usort(XList);
roots([A, B, C, D]) ->
    Ar = reciprocal(A), % Ar = 1/A
    roots([1, mltp(B, Ar), mltp(C, Ar), mltp(D, Ar)]).

roots_filter(Params) ->
    fun(X) ->
        FoldFun = fun(Param, {XpN, Ans}) ->
            {mltp(XpN, X), add(Ans, mltp(Param, XpN))}
        end,
        {_, Complex} = lists:foldr(FoldFun, {1, 0}, Params),
        absolute(Complex) < ?SMALL
    end.
