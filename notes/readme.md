Graph a show original expression E containing:
    - Constants (-1, 1, 2, 4)
    - Symbols (M)
    - BaseScalarField (r)

Graph b shows the result of `E.subs(Symbol('M'), Symbol('r'))`

We can see that it replaced symbol M with r as excpected

Graph c shows the result of `E.subs(Symbol('r'), Symbol('M'))`

We can see that it actually replaced the coordinates of CoordSystem in BaseScalarField
