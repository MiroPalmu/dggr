import dggr
import sympy
from sympy.diffgeom import metric_to_Christoffel_2nd

if __name__=="__main__":
    sympy.init_printing()
    print("Get spacetime patch using: dggr.get_spacetime()")
    print(type(dggr.get_spacetime_patch()))

    print()
    print("Get spacetime coordinate system using get_spacetime_coordsystem(x0, x1, x2, x3)")
    t, x, y, z = sympy.symbols("t x y z")
    coords = dggr.get_spacetime_coordsystem(t, x, y, z)
    print(type(coords))
    print(f"Variables: {coords.symbols}")

    print()
    print("Get Schawrzschild metric in different coordinates: get_schawrzschild_metric(coordinates)")
    g, (t, r, h, phi) = dggr.get_schwarzschild_metric()
    sympy.pprint(t)
    sympy.pprint(r)
    sympy.pprint(h)
    sympy.pprint(phi)
    sympy.pprint(g)
    sympy.pprint(metric_to_Christoffel_2nd(g))

