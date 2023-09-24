import sympy
from sympy import diffgeom as dg

def get_spacetime_patch():
    M = dg.Manifold("Spacetime", 4)
    P = dg.Patch('P', M)
    return P

def get_spacetime_coordsystem(x0, x1, x2, x3):
    return dg.CoordSystem(
            'Spacetime_coords',
            get_spacetime_patch(),
            (x0, x1, x2, x3)
        )
                            
def get_schwarzschild_metric(coordinates: str = "schwarzschild"):
    TP = dg.TensorProduct
    if coordinates == "schwarzschild":
        t_s, r_s, theta_s, phi_s, M = sympy.symbols("t, r, theta, phi, M") 
        coord_system = get_spacetime_coordsystem(t_s, r_s, theta_s, phi_s)
        dt, dr, dtheta, dphi = coord_system.base_oneforms()

        t, r, theta, phi = coord_system.coord_functions()
        g = -(1 - 2 * M / r) * TP(dt, dt) + (1 / (1 - 2 * M / r)) * TP(dr, dr) + r**2 * (TP(dtheta, dtheta) + sympy.sin(theta) ** 2 * TP(dphi, dphi))

        return (g, (t, r, theta, phi)) 
    else if coordinates == "spherical_isotropic":
        t_s, r_s, theta_s, phi_s, M = sympy.symbols("t, r, theta, phi, M") 
        coord_system = get_spacetime_coordsystem(t_s, r_s, theta_s, phi_s)
        dt, dr, dtheta, dphi = coord_system.base_oneforms()

        t, r, theta, phi = coord_system.coord_functions()
        g00 = -((1 - M / (2*r) ) / (1 + M / (2 * r)))**2
        g11 = (1 + M / (2 * r))**4
        g = g00 * TP(dt, dt) + g11 * TP(dr, dr) + r**2 * (TP(dtheta, dtheta) + sympy.sin(theta) ** 2 * TP(dphi, dphi))

        return (g, (t, r, theta, phi)) 

    else:
        raise Exception("Unregonized coordinates for get_schwarzschild_metric!")
from sympy import diffgeom as dg

