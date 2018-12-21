from sympy import *
import numpy as np
import time

# Generates g, first derivatives of g and Christoffel symbols
class g_ginv_Gamma:

    D = None

    x = None
    v = None

    g = None
    dg = None
    Gamma = None

    gf = None
    dgf = None
    Gammaf = None

#    tdot = None
#    U = None

    #=========================
    #=========================
    def __init__( self, x, g ):
        D = len(x)
        v = symbols('v^1:%d'%D)

        print
        print 'Geodesics generator initiated!'
        print '------------------------------'
        # --------------------------------------------------#
        # CONSTRUCT DxDxD matrix OF 1ST DERIVATIVES OF METRIC
        print 'Generating 1st derivaties of g'
        t1 = time.time()
        dg = []
        for rho in range(D):
            dgaux = zeros(D,D)
            for mu in range(D):
                dgaux[mu,mu] = diff(g[mu,mu],x[rho])
                for nu in range(mu+1,D):
                    dgaux[mu,nu] = diff(g[mu,nu],x[rho])
                    dgaux[nu,mu] = dgaux[mu,nu]
            dg.append(dgaux)
#            dg.append(simplify(dgaux))
        t2 = time.time()
        print 'DONE! Time taken:', t2-t1
        # --------------------------------------------------#

        print

        # --------------------------------------------------#
        # CONSTRUCT D DxD MATRICES OF CHRISTOFFEL SYMBOLS
        print 'Generating Christoffel symbols'
        t1 = time.time()
        Christoffel = []
        ginv = g.inv()
#        ginv = simplify(g.inv())
        t3 = time.time()
        print 'Time for g^-1:', t3-t1
        for rho in range(D):
            Gamma = zeros(D,D)
            for mu in range(D):
                for sigma in range(D):
                    Gamma[mu,mu] += ginv[rho,sigma]*(2*diff(g[sigma,mu],x[mu]) - diff(g[mu,mu],x[sigma]))
                for nu in range(mu+1,D):
                    for sigma in range(D):
                        Gamma[mu,nu] += ginv[rho,sigma]*(diff(g[sigma,nu],x[mu]) + diff(g[mu,sigma],x[nu]) - diff(g[mu,nu],x[sigma]))
                    Gamma[nu,mu] = Gamma[mu,nu]
            Christoffel.append(Gamma/2)
#            Christoffel.append(simplify(Gamma/2))
        t2 = time.time()
        print 'DONE! Time taken:', t2-t1
        print
        # --------------------------------------------------#


        # --------------------------------------------------#
        # CONSTRUCT tdot = dt/dtau
#        tdot = g[0,0]
#        for i in range(1,D):
#            tdot += 2*g[0,i] * v[i-1] + g[i,i] * v[i-1]**2
#            for j in range(i+1,D):
#                tdot += 2*g[i,j] * v[i-1]*v[j-1]
#        tdot = simplify(1/sqrt(-tdot))
        # --------------------------------------------------#


        # --------------------------------------------------#
        # CONSTRUCT COVARIANT 4-VELOCITY = tdot*g([1,vec{v}],[1,vec{v}])
#        U = zeros(D)
#        for mu in range(D):
#            U[mu] += g[mu,0]
#            for i in range(1,D):
#                U[mu] += g[mu,i] * v[i-1]
#        U *= tdot
        # --------------------------------------------------#

        self.D = D
        self.x = x
        self.v = v

        self.g = g
        self.dg = dg
        self.Gamma = Christoffel

        self.gf = lambdify(x,g)
        self.dgf = lambdify(x,dg)
        self.Gammaf = lambdify(x,Christoffel)

#        self.tdot = tdot
#        self.U = U
