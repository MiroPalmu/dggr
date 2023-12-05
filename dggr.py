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
    elif coordinates == "spherical_isotropic":
        t_s, r_s, theta_s, phi_s, M = sympy.symbols("t, r, theta, phi, M") 
        coord_system = get_spacetime_coordsystem(t_s, r_s, theta_s, phi_s)
        dt, dr, dtheta, dphi = coord_system.base_oneforms()

        t, r, theta, phi = coord_system.coord_functions()
        g00 = -((1 - M / (2*r) ) / (1 + M / (2 * r)))**2
        g11 = (1 + M / (2 * r))**4
        g = g00 * TP(dt, dt) + g11 * TP(dr, dr) + r**2 * (TP(dtheta, dtheta) + sympy.sin(theta) ** 2 * TP(dphi, dphi))

        return (g, (t, r, theta, phi)) 

    elif coordinates == "kerr-schild":
        x0_s, x1_s, x2_s, x3_s, M = sympy.symbols("x0, x1, x2, x3, M")
        coord_system = get_spacetime_coordsystem(x0_s, x1_s, x2_s, x3_s)
        dx0, dx1, dx2, dx3 = coord_system.base_oneforms()

        x0, x1, x2, x3 = coord_system.coord_functions()
        r = sympy.sqrt(x0**2 + x1**2 + x2**2 + x3**2)
        H = M / r

        l_oneform = dx0 + (x1 / r) * dx1 + (x2 / r) * dx2 + (x3 / r) * dx3
        minkowski_metric = -TP(dx0, dx0) + TP(dx1, dx1) + TP(dx2, dx2) +  TP(dx3, dx3)

        g = minkowski_metric + 2 * H * TP(l_oneform, l_oneform)

        return (g, (x0, x1, x2, x3))

    else:
        raise Exception("Unregonized coordinates for get_schwarzschild_metric!")


def get_euclidean_patch(n: int):
    M = dg.Manifold("Spacetime", n)
    P = dg.Patch('P', M)
    return P

def get_euclidean_coordsystem(*coordinates):
    return dg.CoordSystem(
            'euclidean_coords',
            get_euclidean_patch(len(coordinates)),
            coordinates
        )
                            
def get_euclidean_metric(coordinates: str = "orthogonal"):
    TP = dg.TensorProduct
    if coordinates == "orthogonal":
        x_s, y_s, z_s = sympy.symbols("x, y, z") 
        coord_system = get_euclidean_coordsystem(x_s, y_s, z_s)
        dx, dy, dz = coord_system.base_oneforms()

        x, y, z = coord_system.coord_functions()
        g = TP(dx, dx) + TP(dy, dy) + TP(dz, dz)

        return (g, (x, y, z)) 

    elif coordinates == "spherical":
        r_s, theta_s, phi_s = sympy.symbols("r, theta, phi") 
        coord_system = get_euclidean_coordsystem(r_s, theta_s, phi_s)
        dr, dtheta, dphi = coord_system.base_oneforms()

        r, theta, phi = coord_system.coord_functions()
        g11 = 1 
        g22 = r**2
        g33 = r**2 * sympy.sin(theta)**2
        g = g11 * TP(dr, dr) + g22 * TP(dtheta, dtheta) + g33 * TP(dphi, dphi)
        return (g, (r, theta, phi)) 

    else:
        raise Exception("Unregonized coordinates for get_euclidean_metric!")
