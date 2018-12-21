from __init__ import *

def metric( D, M = Symbol('M') ):
    dSph = DSphere(D-1)

    t,r,pi = symbols('t,r,pi')

    g = Matrix(D,D, lambda i,j: 0)
    f1 = simplify(8*pi*M/((D-2)*dSph.AreaDm1(r)))
    f = 1-2*r*f1
    g[0,0] = -f
    g[1,1] = 1/f
    gSphere,theta = dSph.gThetaTheta()
    for i in range(2,D):
        g[i,i] = r**2*gSphere[i-2]

    x=[]
    x.append(t)
    x.append(r)
    for i in range(D-2):
        x.append(theta[i])

    return g.subs({pi:np.pi}),x
#======================================================
def eventHorizon(D,M):
    dSph = DSphere(D-1)
    RS = 16*pi*M/((D-2)*dSph.AreaDm1(1.0))
    rH = RS.subs({pi:np.pi})
    return rH**(1.0/(D-3))
#======================================================
