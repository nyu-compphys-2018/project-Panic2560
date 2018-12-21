from sympy import *
import numpy as np
import time

# Geodesics solver. Takes functions of metric, Christoffel symbols and extra
# forces and then, after determining initial conditions, outputs evolution of
# coordinate position, coordinate velocity and D-velocity.

class GEODESICS_SOLVER:

    g = None
    dg = None
    Gamma = None

    XFORCE = None
    XFargs = None

    t0 = None
    x0 = None
    v0 = None

    t = None
    x = None
    v = None

    N = None

    INTEGRATOR = None
    INTEGRATOR_ERROR_ORDER = None
    #=========================
    #=========================
    def __init__( self, gf, dgf, Gammaf, XFORCE ):

        self.g = gf
        self.dg = dgf
        self.Gamma = Gammaf
        self.XFORCE = XFORCE

    #=========================
    #=========================
    def SetInitialConditions( self, t0, x0, v0 ):
        self.t0 = float(t0)
        self.x0 = np.array( x0 )
        self.v0 = np.array( v0 )
    #=========================
    def SetXForceArgs( self, *XFargs ):
        self.XFargs = XFargs
    #=========================
    def SetIntegrator( self, INTEGRATOR, INTEGRATOR_ERROR_ORDER ):
        self.INTEGRATOR = INTEGRATOR
        self.INTEGRATOR_ERROR_ORDER = INTEGRATOR_ERROR_ORDER
    #=========================
    #=========================
    def FORCE( self, t, x, v ):
        d = len(x)

        f = np.zeros(d, float)
        y = np.zeros(d+1, float)


        y[0] = t
        y[1:] = x

        Xftmp = self.XFORCE( t, x, v, *self.XFargs )

        for i in range(d):
            f[i] += self.Gamma(*y)[0][0][0] * v[i] - self.Gamma(*y)[i+1][0][0]
            for j in range(d):
                f[i] += 2*(self.Gamma(*y)[0][j+1][0] * v[j]*v[i] -  self.Gamma(*y)[i+1][j+1][0] * v[j])
                for k in range(d):
                    f[i] += self.Gamma(*y)[0][j+1][k+1] * v[j]*v[k]*v[i] - self.Gamma(*y)[i+1][j+1][k+1] * v[j]*v[k]
            f[i] += ( Xftmp[i+1] - Xftmp[0]*v[i] )/( self.DtDtau( t, x, v ) )**2

        return f
    #=========================
    def evolve( self, tN, N, SYMPLECTIC_FLAG = None ):

        self.N = N
        if SYMPLECTIC_FLAG == True:  
            self.t,self.x,self.v = self.INTEGRATOR( self.t0, self.x0, self.v0, tN, N, self.FORCE, *[] )
        else:
            d = len(self.x0)

            dt = (tN - self.t0)/N

            t = self.t0 + dt*np.arange(N)
            x = np.zeros([d,N], float)
            v = np.zeros([d,N], float)

            x[:,0] = self.x0
            v[:,0] = self.v0

            for n in range(N-1):
                x[:,n+1],v[:,n+1] = self.INTEGRATOR( t[n], x[:,n], v[:,n], dt, self.FORCE, *[] )

            self.t = t
            self.x = x
            self.v = v
    # =============================================#
    def evolveAdaptive( self, tN, sigma, dt0 ):
	d = len(self.x0)

        t = []
        x = []
	v = []

        taux = self.t0
	xaux = self.x0
        vaux = self.v0

	t.append(self.t0)
        x.append(self.x0)
        v.append(self.v0)

        l = self.INTEGRATOR_ERROR_ORDER
        lM1 = l-1
        lFactor = 2**l - 2

        dt = dt0
        N = 1
        while taux < tN:
            N += 1
            x1aux,v1aux = self.INTEGRATOR( taux, xaux, vaux, dt, self.FORCE, *[] )
            x1,v1 = self.INTEGRATOR( taux+dt, x1aux, v1aux, dt, self.FORCE, *[] )
            x2,v2 = self.INTEGRATOR( taux, xaux, vaux, 2*dt, self.FORCE, *[] )
            Dxv2 = sum((x1-x2)**2) + sum((v1-v2)**2)
            rho = sigma*dt*lFactor/np.sqrt(Dxv2)
            if rho <= 1.0:
                taux += dt
                xaux,vaux = x1aux,v1aux
                dt = dt*rho**(1.0/lM1)
            else:
                dt = dt*rho**(1.0/lM1)
                xaux,vaux = self.INTEGRATOR( taux, xaux, vaux, dt, self.FORCE, *[] )
                taux += dt
            t.append(taux)
            x.append(xaux)
            v.append(vaux)

        xArr = np.zeros([d,N], float)
        vArr = np.zeros([d,N], float)
        for i in range(d):
            for n in range(N):
                xArr[i,n] = x[n][i]
                vArr[i,n] = v[n][i]

        self.t = np.array(t, float)
        self.x = xArr
        self.v = vArr
        self.N = N
    #=========================
    def DtDtau(self, t, x, v):
        D = len(x) + 1

        y = np.zeros(D, float)

        y[0] = t
        y[1:] = x

        tdot = self.g(*y)[0][0]
        for i in range(1,D):
            tdot += 2*self.g(*y)[0][i] * v[i-1] + self.g(*y)[i][i] * v[i-1]**2
            for j in range(i+1,D):
                tdot += 2*self.g(*y)[i][j] * v[i-1]*v[j-1]
        return 1/np.sqrt(-tdot)
    #=========================
    def FourMomenta( self ):
        D = len(self.x0)+1

        U = np.zeros([D,self.N], float)
        y = np.zeros(D, float)

        for n in range(self.N):
            y[0] = self.t[n]
            y[1:] = self.x[:,n]

            for mu in range(D):
                U[mu,n] += self.g(*y)[mu][0]
                for i in range(1,D):
                    U[mu,n] += self.g(*y)[mu][i] * self.v[i-1,n]

            U[:,n] =  U[:,n] * self.DtDtau( self.t[n], self.x[:,n], self.v[:,n] )
        return U
    #=========================
    def v0i_umuConstraint( self, i, mu, umu ):
        # Gives vi given vj (j not equal to i) and and the mu-th component of the covariant 4-velocity
        # Requirement: g_nu,i are non-zero for at least one value of nu
        d = len(self.x0)
        y = np.zeros(d+1, float)
        y[0] = self.t0
        y[1:] = self.x0
        g = self.g(*y)

        #==========================
        a = g[i+1][i+1]*umu**2 + (g[mu][i+1])**2

        b = g[i+1][0]*umu**2 + g[mu][i+1]*g[mu][0]
        for j in range(i) + range(i+1,d):
            b += ( g[i+1][j+1]*umu**2 + g[mu][i+1]*g[mu][j+1] ) * self.v0[j]
        b *= 2

        c1 = g[0][0]
        c2 = g[mu][0]
        for j in range(i) + range(i+1,d):
            c1 += 2*g[0][j+1] * self.v0[j]
            c2 += g[mu][j+1] * self.v0[j]
            for k in range(i) + range(i+1,d):
                c1 += g[j+1][k+1] * self.v0[j]*self.v0[k]
        c = c1*umu**2 + c2**2

        #==========================
        tdot0 = g[0][0]
        for j in range(i) + range(i+1,d):
            tdot0 += 2*g[0][j+1] * self.v0[j]
            for k in range(i) + range(i+1,d):
                tdot0 += g[j+1][k+1] * self.v0[j]*self.v0[k]

        tdot1 = g[0][i+1]
        for j in range(i) + range(i+1,d):
            tdot1 += g[i+1][j+1] * self.v0[j]
        tdot1 *= 2

        tdot2 = g[i+1][i+1]
        #==========================
        Umu = g[mu][0]
        for j in range(i) + range(i+1,d):
            Umu += g[mu][j+1] * self.v0[j]
        #==========================
        Discr = b**2 - 4*a*c
        if umu < 1.e-10:
            vi = -b/2/a
            mtdot2 = tdot0 + tdot1*vi + tdot2*vi**2
            Umu = (Umu + g[mu][i+1]*vi)/np.sqrt(-mtdot2)
            print '+/-:', vi, mtdot2
            self.v0[i] = vi
            print 'U[',mu,']_+ = ', Umu
        elif Discr<0:
            print 'No real roots'
        else:
            vip = (-b+np.sqrt(Discr))/2/a
            vim = (-b-np.sqrt(Discr))/2/a

            mtdotp2 = tdot0 + tdot1*vip + tdot2*vip**2
            mtdotm2 = tdot0 + tdot1*vim + tdot2*vim**2

            Umup = (Umu + g[mu][i+1]*vip)/np.sqrt(-mtdotp2)
            Umum = (Umu + g[mu][i+1]*vim)/np.sqrt(-mtdotm2)

            print '+:', vip, mtdotp2
            print '-:', vim, mtdotm2
            print umu, Umup, Umum
            if abs(umu-Umup) < 1.e-5:
                self.v0[i] = vip
                print 'U[',mu,']_+ = ', Umup
            elif abs(umu-Umum) < 1.e-5:
                self.v0[i] = vim
                print 'U[',mu,']_- = ', Umum
    #==============================================
    def v0All_uConstaint( self, u, TOL=1.e-6 ):
        d = len(self.x0)
        D = d+1
        y = np.zeros(D, float)
        y[0] = self.t0
        y[1:] = self.x0
        g = self.g(*y)

        #==========================
        v = symbols('v1:%d'%D)
        ftmp = zeros(D)
        f = zeros(d)

        for mu in range(D):
            ftmp[mu] += g[0][0]*u[mu]**2 + g[mu][0]**2
            for i in range(1,D):
                ftmp[mu] += 2*( g[0][i]*u[mu]**2 + g[mu][0]*g[mu][i] ) * v[i-1]
                for j in range(1,D):
                    ftmp[mu] += ( g[i][j]*u[mu]**2 + g[mu][i]*g[mu][j] ) * v[i-1]*v[j-1]

        f[0] = ftmp[0]
        for i in range(2,D): # Exclude u1
            f[i-1] = simplify(ftmp[i])

        J = Matrix(d,d, lambda i,j: diff(f[i],v[j]))

        t1 = time.time()
        Jinv = simplify(J.inv())
        t2 = time.time()
        print t2-t1
        #==========================
        ff = lambdify(v,f)
        Jinvf = lambdify(v,Jinv)

        v0 = 1.e-10*np.ones(d, float)
        v0tmp = v0 + 2*TOL
        Niter = 0
        while abs(v0[0]-v0tmp[0])>TOL:
            Niter += 1
            v0tmp = v0
            for i in range(d):
                factor = 0.0
                for j in range(d):
                    factor += Jinvf(*v0)[i][j]*ff(*v0)[0][j]
                v0[i] -= factor

        print 'Done after', Niter, 'iterations'
        print 'E/m = ', ff(*v0)[0][0]
        print 'L/m = ', ff(*v0)[0][d-1]

        #==============================================
#        v0[0] = 1.e-8
        self.v0 = v0
    #=========================
#============================================================================
#============================================================================
#============================================================================
#============================================================================
#============================================================================
# These here could be useful for later development. They are symbolic
# coordinate transformations and a root finder for the event horizon of
# black holes
#----------------------------------------------------------------------
# Type 1: xold = f(xnew)
# Type 2: xnew = f(xold)
#def DiffTransf( xold, xnew, typeOpt ):
#    D = len(xold)
#    if typeOpt == 1:
#        Diff1Form = Matrix(D,D,lambda mu,nu: diff(xold[mu],xnew[nu]))
#        DiffVector = simplify(Diff1Form.inv())
#    elif typeOpt == 2:
#        DiffVector = Matrix(D,D,lambda mu,nu: diff(xnew[mu],xold[nu]))
#        Diff1Form = simplify(DiffVector.inv())
#    return DiffVector, Diff1Form

#=========================
#def eventHorizon( self, r0=3.0, TOL=1.e-12 ):
#    # WORKS ONLY FOR r-dependent g_11 = g_11(r) and g_mu1 = g_11 delta_1^mu
#    g11 = lambdify(self.x[1],1/self.g[1,1])
#    dg11 = lambdify(self.x[1],diff(1/self.g[1,1],self.x[1]))
#    r = r0 + 1.e-5
#    rprime = r0
#    N = 0
#    while abs(r-rprime)>TOL and N<1000:
#        N += 1
#        r = rprime
#        rprime = r - g11(r)/dg11(r)
#
#    if N == 10000:
#        print 'NO EVENT HORIZON FOUND AFTER 10000 ITERATIONS'
#        return rprime
#    else:
#        return rprime
