-module(linalg_roots).
-vsn('1.0').
-author('simon.klassen').
-import(math,[pi/0,acos/1,cos/1,sqrt/1,pow/2]).
-export([roots/1]).
-define(SMALL,1.0e-10).

-type vector() :: list(float).
-spec roots(vector()) -> vector().
% null case
roots([])->
    [];

% singleton
roots([_])->
    [];

% ax + b = 0
roots([A,B])->
    [-B/A];

% ignore higher order when alpha is small
roots([A|Tail]) when abs(A)<?SMALL ->
    roots(Tail);

% quadratic polynomial
% well known closed form
% ax^2 + bx + c = 0 
roots([A,B,C]) ->
    Q = -0.5*(B + sign(1, B) * sqrt(B*B - 4*A*C)),
    X1 = Q/A,
    X2 = C/Q,
    [min(X1,X2),max(X1,X2)];

% cubic polynomial
% viete
% x^3 + ax^2 + bx +c = 0
roots([1,A,B,C]) ->
    Q = (A*A-3*B)/9,
    R = (2*A*A*A-9*A*B+27*C)/54,

    Q3 = Q*Q*Q,
    R2 = R*R,

    case R2 < Q3 of
      true-> 
        % three real roots
        TH=acos(R / sqrt(Q3)),
        SqrtQ=sqrt(Q),
        X1=-2*SqrtQ*cos(TH/3) - A/3,
        X2=-2*SqrtQ*cos((TH + 2*pi())/3) - A/3,
        X3=-2*SqrtQ*cos((TH - 2*pi())/3) - A/3,
        lists:sort([X1,X2,X3]);
      false-> 
        % one real root
        Alpha = -sign(1,R)*pow(abs(R)+sqrt(R2-Q3),1/3),
        Beta = case Alpha of
           Alpha when abs(Alpha)<?SMALL -> 0;
           Alpha -> Q/Alpha
         end,
        [(Alpha + Beta) - A/3]
    end.

sign(A,B)->
    case B>=0.0 of
        true->abs(A);
        false-> -abs(A)
    end.
