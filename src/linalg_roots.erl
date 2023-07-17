-module(linalg_roots).
-vsn('1.0').
-author('simon.klassen').
-import(math, [pi/0, acos/1, cos/1, pow/2]).
-import(linalg_complex, [sqrt/1, qbrt/1]).
-import(linalg_complex, [add/1, add/2, mltp/1, mltp/2, reciprocal/1]).
-export([roots/1]).
-define(SMALL, 1.0e-10).

-type vector() :: list(number()).
-type complex() :: Real::number() | {Real::number(), Imaginary::number()}.
-spec roots(vector()) -> complex().
% null case
roots([]) ->
    [];
% singleton
roots([_]) ->
    [];
% ignore higher order when alpha is small
roots([A | Tail]) when abs(A) < ?SMALL ->
    roots(Tail);
% ax + b = 0
roots([A, B]) ->
    [-B / A];
% quadratic polynomial
% well known closed form
% ax^2 + bx + c = 0
roots([A, B, C]) ->
    Denom = reciprocal(mltp(2,A)), % 1/(2*A)
    P = mltp(Denom, -B), % -b/(2*A)
    % Q = sqrt(B^2-4*A*C)/(2*A)
    Q = mltp(Denom, sqrt(add(mltp(B, B), mltp([-4, A, C])))),
    lists:usort([add(P, Q), add(P, mltp(-1, Q))]); % P +/- Q
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
    
    % Wp = (-1+sqrt(3)i)/2
    % Wn = (-1-sqrt(3)i)/2
    Wp = {-1/2, sqrt(3)/2},
    Wn = {-1/2, -sqrt(3)/2},
    
    % QpR = qbrt(Q+R)
    % QnR = qbrt(Q-R)
    QpR = qbrt(add(Q, R)),
    QnR = qbrt(add(Q, mltp(-1, R))),
    
    % X1 = P + QpR + QnR
    % X2 = P + Wp * QpR + Wn * QnR
    % X3 = P + Wn * QpR + Wp * QnR
    X1 = add([P, QpR, QnR]),
    X2 = add([P, mltp(Wp, QpR), mltp(Wn, QnR)]),
    X3 = add([P, mltp(Wn, QpR), mltp(Wp, QnR)]),
    lists:usort([X1, X2, X3]);
roots([A, B, C, D]) ->
  roots([1, B/A, C/A, D/A]).

sign(A, B) ->
    case B >= 0.0 of
        true -> abs(A);
        false -> -abs(A)
    end.
