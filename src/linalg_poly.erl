-module(linalg_poly).
-export([polyval/2,polyfit/3]).

polyval([A],_X)->
    A;
polyval([A,B],X)->
    A + B * X;
polyval([A,B,C],X)->
    A + B * X + C * X * X.

polyfit(_Xs,Ys,0)->
    [linalg:mean(Ys)];

polyfit(Xs,Ys,1)->
    XYs = lists:zip(Xs,Ys),
    N = length(XYs),

    {Xhat, Yhat} = lists:foldl(fun({X,Y}, {S1,S2}) -> {S1+X/N,S2+Y/N} end, {0,0}, XYs),

    {Sxx, Sxy} = lists:foldl(fun({X,Y}, {S1,S2}) -> 
                                     {S1+(X-Xhat)*(X-Xhat),
                                      S2+(X-Xhat)*(Y-Yhat)} 
                             end, {0,0}, XYs),

    B = Sxy / Sxx,
    A = Yhat - B * Xhat,

    [A,B];

polyfit(Xs,Ys,2)->
    XYs = lists:zip(Xs,Ys),
    N = length(XYs),

    {Xhat, X2hat, Yhat} = lists:foldl(fun({X,Y}, {S1,S2,S3}) -> {S1+X/N,S2+X*X/N,S3+Y/N} end, {0,0,0}, XYs),

    {Sxx, Sxy, Sxx2, Sx2x2, Sx2y} = lists:foldl(fun({X,Y}, {S1,S2,S3,S4,S5}) -> 
                                                        {S1+(X-Xhat)*(X-Xhat),
                                                         S2+(X-Xhat)*(Y-Yhat),
                                                         S3+(X-Xhat)*(X*X-X2hat),
                                                         S4+(X*X-X2hat)*(X*X-X2hat),
                                                         S5+(X*X-X2hat)*(Y-Yhat)} 
                                                end, {0,0,0,0,0}, XYs),

    B = (Sxy * Sx2x2 - Sx2y * Sxx2) / (Sxx * Sx2x2 - Sxx2 * Sxx2),
    C = (Sx2y * Sxx - Sxy * Sxx2) / (Sxx * Sx2x2 - Sxx2 * Sxx2),
    A = Yhat - B * Xhat - C * X2hat,

    [A,B,C].
