from sympy import *
import numpy as np
import time

# Produces initial conditions for particle 2 and masses of particles 2 and 3
# after the decay of parent particle given vf1 and v03.

# Running this python script directly solve for the final constants of motion
# without having to solve the geodesics using PenroseProcess.py

class PENROSE_SPLIT_IC:

    g = None

    m1 = None
    m2 = None
    m3 = None

    v1 = None
    v2 = None
    v3 = None

    u1contra = None
    u2contra = None
    u3contra = None

    u1cov = None
    u2cov = None
    u3cov = None

    #=================================
    #=================================
    def __init__( self, g, m1, v1 ):
        self.g = g
        self.m1 = m1
        self.v1 = v1
        self.u1contra = self.ucontra( v1 )
        self.u1cov = self.ucov( v1 )
        print
        print 'v1 = ', self.v1
        print 'u_1mu = ', self.u1cov
        print
    #=================================
    def Setv2( self, v2 ):
        self.v2 = v2
        self.u2contra = self.ucontra( v2 )
        self.u2cov = self.ucov( v2 )
        print
        print 'v2 = ', self.v2
        print 'u_2mu = ', self.u2cov
        print
    def Setv3( self, v3 ):
        self.v3 = v3
        self.u3contra = self.ucontra( v3 )
        self.u3cov = self.ucov( v3 )
        print
        print 'v3 = ', self.v3
        print 'u_3mu = ', self.u3cov
        print
    #=================================
    def m2m3v2d( self ):
        d = len(self.v1)
        D = d+1

        u1u3 = 0.0
        for mu in range(D):
            for nu in range(D):
                u1u3 += self.g[mu][nu] * self.u1cov[mu]*self.u3cov[nu]

        m3s = Symbol('m3')
        m2s = sqrt(1 + m3s**2 + 2*m3s*u1u3)

        L2 = (self.u1cov[d] - m3s* self.u3cov[d])/m2s
        v2s = zeros(1,d)

        t1dot = self.DtDtau(self.v1)
        lambda12 = (t1dot*self.v1[0])**2
        k = L2**2*(1 - lambda12/(m2s**2 + lambda12*self.g[1][1]))/self.g[0][0]
        psi = self.g[0][0]*self.g[d][d] - (self.g[0][d])**2

#        t1 = time.time()
        v2s[d-1] = -self.g[0][d]*(1+k) - sign(self.u3cov[d])*sqrt(-k*psi*(1+k))
        v2s[d-1] = simplify(self.g[0][0]*v2s[d-1]/((self.g[0][d])**2 + k*self.g[0][0]*self.g[d][d]))
#        t2 = time.time()
#        print 'Time for v2d(m3): Dt = ', t2-t1

        v2s[0] = simplify(sqrt(-(self.g[0][0] + 2*self.g[0][d]*v2s[d-1] + self.g[d][d]*v2s[d-1]**2)*lambda12/(m2s**2 + lambda12*self.g[1][1])))
#        t1 = time.time()
#        print 'Time for v20(m3): Dt = ', t1-t2

        for i in range(1,d-1):
            v2s[i] = 0

        E1 = -self.u1cov[0]
        E3 = -self.u3cov[0]
#        print E1,E3

        E2 = simplify((E1 - m3s*E3)/m2s)

        # Defin3 f0 and use NR to find m3 and, thus, m2 and v2d
        f1 = self.g[0][0]
        f2 = self.g[0][0]
        for i in range(d):
            f1 += self.g[0][i+1]*v2s[i]
            f2 += 2*self.g[0][i+1]*v2s[i]
            for j in range(d):
                f2 += self.g[i+1][j+1]*v2s[i]*v2s[j]
        f0 = f1**2 + f2*E2**2

        df0 = diff(f0,m3s)
#        t1 = time.time()
        ftmp = f0/df0
#        t2 = time.time()
#        print 'Time to construct f(m3)=0: Dt = ', t2-t1
#        print 'f(m3) = ', ftmp

        f = lambdify(m3s,ftmp)
        m2 = lambdify(m3s,m2s)
        v2 = lambdify(m3s,v2s)

        m3 = -1.e-8
        m3tmp = m3 + 2.e-8
        Niter = 0
        while abs(m3-m3tmp)>1.e-8 and Niter<10000:
            Niter += 1
            m3tmp = m3
            m3 -= f(m3)
#            print Niter, m3, m3tmp, abs(m3-m3tmp)

        v2self = v2(m3)[0]
        self.m2 = self.m1*m2(m3)
        self.m3 = self.m1*m3
        self.v2 = v2self
        self.u2contra = self.ucontra( v2self )
        self.u2cov = self.ucov( v2self )

        print
        print 'm: ', self.m1, self.m2, self.m3
        print 'v: ', self.v1, self.v2, self.v3
        print 'u_mu: ', self.u1cov, self.u2cov, self.u3cov
        print 'Dm/mmax = ', (self.m1 - self.m2 - self.m3)/max(self.m1,self.m2,self.m3)
        print 'Dp_mu = ', self.m1*self.u1cov - self.m2*self.u2cov - self.m3*self.u3cov
        print

    #=================================
    def DtDtau( self, v ):
        D = len(v) + 1
        tdot = self.g[0][0]
        for i in range(1,D):
            tdot += 2*self.g[0][i] * v[i-1]
            for j in range(1,D):
                tdot += self.g[i][j] * v[i-1]*v[j-1]
        tdot = 1/np.sqrt(-tdot)
        return tdot
    #=================================
    def ucontra( self, v ):
        D = len(v)+1

        u = np.zeros(D, float)

        u[0] = 1.0
        for i in range(1,D):
            u[i] = v[i-1]

        return u*self.DtDtau(v)
    #=================================
    def ucov( self, v ):
        D = len(v)+1

        u = np.zeros(D, float)

        for mu in range(D):
            u[mu] += self.g[mu][0]
            for i in range(1,D):
                u[mu] += self.g[mu][i] * v[i-1]

        return u*self.DtDtau(v)
#===============================================================================
def v0i_umuConstaint( i, mu, umu, g, v ):
    d = len(v)
    #==========================
    a = g[i+1][i+1]*umu**2 + g[mu][i+1]**2

    b = g[0][i+1]*umu**2 + g[mu][i+1]*g[mu][0]
    for j in range(i) + range(i+1,d):
        b += ( g[i+1][j+1]*umu**2 + g[mu][i+1]*g[mu][j+1] ) * v[j]
    b *= 2

    c1 = g[0][0]
    c2 = g[mu][0]
    for j in range(i) + range(i+1,d):
        c1 += 2*g[0][j+1] * v[j]
        c2 += g[mu][j+1] * v[j]
        for k in range(i) + range(i+1,d):
            c1 += g[j+1][k+1]*umu**2 * v[j]*v[k]
    c = c1*umu**2 + c2**2
    #==========================
    tdot0 = g[0][0]
    for j in range(i) + range(i+1,d):
        tdot0 += 2*g[0][j+1] * v[j]
        for k in range(i) + range(i+1,d):
            tdot0 += g[j+1][k+1] * v[j]*v[k]

    tdot1 = g[0][i+1]
    for j in range(i) + range(i+1,d):
        tdot1 += g[i+1][j+1] * v[j]
    tdot1 *= 2

    tdot2 = g[i+1][i+1]
    #==========================
    Umu = g[mu][0]
    for j in range(i) + range(i+1,d):
        Umu += g[mu][j+1] * v[j]
    #==========================
    Discr = b**2 - 4*a*c
    if abs(umu) < 1.e-10:
        vi = -b/2/a
        mtdot2 = tdot0 + tdot1*vi + tdot2*vi**2
        Umu = (Umu + g[mu][i+1]*vi)/np.sqrt(-mtdot2)
        print
        print '+/-:', vi, mtdot2
        print 'U[',mu,']_+/- = ', Umu
        print
        return vi
    elif Discr<0:
        print
        print 'No real roots'
        print
        return 0.0
    else:
        vip = (-b+np.sqrt(Discr))/2/a
        vim = (-b-np.sqrt(Discr))/2/a

        mtdotp2 = tdot0 + tdot1*vip + tdot2*vip**2
        mtdotm2 = tdot0 + tdot1*vim + tdot2*vim**2

        Umup = (Umu + g[mu][i+1]*vip)/np.sqrt(-mtdotp2)
        Umum = (Umu + g[mu][i+1]*vim)/np.sqrt(-mtdotm2)

        print
        print '+:', vip, mtdotp2
        print '-:', vim, mtdotm2
        if abs(umu-Umup) < 1.e-5:
            print 'U[',mu,']_+ = ', Umup
            print
            return vip
        elif abs(umu-Umum) < 1.e-5:
	    print 'U[',mu,']_- = ', Umum
            print
            return vim
    #==============================================            
if __name__=='__main__':

    d = 3

    g = np.zeros([d+1,d+1], float)
    g[0][0] = 0.1612838
    g[0][3] = g[3][0] = -1.04515542
    g[1][1] = 8.94425095
    g[2][2] = 2.96608268
    g[3][3] = 4.71672256

    print 'g(x) = '
    print g
    print
    v1 = np.zeros(d, float)
    v2 = np.zeros(d, float)
    v3 = np.zeros(d, float)

    m1 = 1.0
    v1[0] = -0.0854537888
    v1[1] = -0.000000
    v1[2] = 0.221585096

    v3[0] = 0.0
    v3[1] = 0.0
    omega = -g[0][0]/g[0][d]
#    psi = np.sqrt(-g[0][0]*g[d][d]+(g[0][d])**2)/g[d][d]
#    print -psi, +psi
#    v3[2] = omega - 1.e-6
    v3[2] = omega - 5.4820298e-2

    q = PENROSE_SPLIT_IC( g, m1, v1 )
    q.Setv3( v3 )

    q.m2m3v2d( )
