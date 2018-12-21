from __init__ import *
#======================================================
def metric( M = Symbol('M'), J = Symbol('J') ):
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
    Delta = r**2 + a**2 - Rs*r

    g[0,0] = -1 + Rs*r/Sigma
    g[1,1] = Sigma/Delta
    g[2,2] = Sigma
    g[3,3] = (r**2 + a**2 + Rs*r*a**2/Sigma * sin(theta)**2)*sin(theta)**2
    g[0,3] = g[3,0] = -Rs*r*a*sin(theta)**2/Sigma

    return g, x
#======================================================
def eventHorizon(M,J):
    a = J/M
    return M + np.sqrt(M**2 - a**2)
#======================================================
def ergoSurface(theta,M,J):
    a = J/M
    return M + np.sqrt(M**2 - a**2*np.cos(theta)**2)
