from __init__ import *
#======================================================
def metric( M = Symbol('M'), J = Symbol('J'), Q = Symbol('Q') ):
    t, r, theta, phi = symbols('t, r, theta, phi')
    x = []
    x.append(t)
    x.append(r)
    x.append(theta)
    x.append(phi)

    g = Matrix(4,4, lambda i,j: 0)
    Rs = 2*M
    a = J/M
    Sigma = r**2 + a**2 * cos(theta)**2
    Delta = r**2 + a**2 - Rs*r + Q**2

    g[0,0] = simplify(a**2*sin(theta)**2 -Delta)/Sigma
    g[1,1] = simplify(Sigma/Delta)
    g[2,2] = simplify(Sigma)
    g[3,3] = simplify(((r**2 + a**2)**2 - Delta*a**2)/Sigma * sin(theta)**2)
    g[0,3] = g[3,0] = (a*Delta*sin(theta)-a*(r**2+a**2)*sin(theta)**2)/Sigma

    return g, x
#======================================================
def eventHorizon( M, J, Q ):
    a = J/M
    return M + np.sqrt(M**2 - a**2 - Q**2)
#======================================================
def ergoSurface( theta, M, J, Q):
    a = J/M
    return M + np.sqrt(M**2 - a**2*np.cos(theta)**2 - Q**2)
