from __init__ import *

def metric( D, M = Symbol('M'), Q = Symbol('Q')):
    dSph = DSphere(D-1)

    t,r,pi = symbols('t,r,pi')

    g = Matrix(D,D, lambda i,j: 0)
    f1 = simplify(-16*pi*M/((D-2)*dSph.AreaDm1(r)))
    f2 = simplify(8*pi*Q**2/((D-2)*(D-3)))
    f = 1 + r*f1 + f2/r**(2*(D-3))
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
def eventHorizon(D,M,Q):
    dSph = DSphere(D-1)
    RS = 16*np.pi*M/((D-2)*dSph.AreaDm1(1.0))
    A = 16*np.pi*Q/((D-2)*(D-3))
    rH = (RS + np.sqrt(RS**2 - A))/2
#    rH1 = (RS + np.sqrt(RS**2 - A))/2
#    rH = rH1.subs({pi:np.pi})
    return rH**(1.0/(D-3))
#======================================================
