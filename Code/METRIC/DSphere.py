import numpy as np
from sympy import *

# Generator of quantities associated with a D-dimensional sphere. Options:
# [1] Spherical-to-Cartesian transformation x_i = r * Omega_i
# [2] Metric g_{ij} on unit sphere (for general radius simply multiply with r^2)
# [3] Volume VolD of D-dimensional compact sphere
# [4] Surface area AreaDm1 of D-dimensional compact sphere
class DSphere:
    d = None
    theta = None

    # =============================================
    def __init__(self, d):
        self.d = d
        self.theta = symbols('theta_2:%d'%(d+1))
    # =============================================
    # =============================================
    def Omega(self, i):
        if i == 1:
            return cos(self.theta[0])
        elif i == self.d:
            omega = 1
            for j in range(self.d-1):
                omega = omega*sin(self.theta[j])
            return omega
        else:
            omega = 1
            for j in range(i-1):
                omega = omega*sin(self.theta[j])
            return omega*cos(self.theta[i-1])
    # =============================================
    def gThetaTheta(self):
        g = []
        for i in range(1,self.d):
            gi = 1
            for j in range(i-1):
                gi = gi*sin(self.theta[j])**2
            g.append(gi)
        return g,self.theta
    # =============================================
    def SphericalCoordinates(self):
        r = Symbol('r')
        y = []
        for i in range(self.d):
            y.append(r*self.Omega(i+1))
        return y
    # =============================================
    def VolD(self,r):
#        pi = np.pi
        pi = Symbol('pi')
        return simplify(r**(self.d)*pi**(self.d/2.0)/Gamma(self.d/2.0+1))/pi*pi
    # =============================================
    def AreaDm1(self,r):
        pi = np.pi
#        pi = Symbol('pi')
        return r**(self.d-1)*self.d*pi**(self.d/2.0)/Gamma(self.d/2.0+1)/pi*pi
#        return simplify(r**(self.d-1)*self.d*pi**(self.d/2.0)/Gamma(self.d/2.0+1))/pi*pi


# =============================================
# =============================================
def Gamma(n):
    if (2*n)%2 != 0:
        pi = np.pi
        Gamma = np.sqrt(np.pi)
#        pi = Symbol('pi')
#        Gamma = sqrt(pi)
        for m in range(1,int(2*n),2):
            Gamma = m*Gamma/2
    else:
        Gamma = 1
        for m in range(1,int(n)):
            Gamma = m*Gamma

    return Gamma
# =============================================
# =============================================
if __name__=='__main__':
    d = int(input('NUmber of spatial dimensions: d = '))
    x = symbols('x1:%d'%(d+1))
    y = DSphere(d).SphericalCoordinates()
    print 'Spherical Coordinates: '
    for i in range(d):
        print x[i], ' = ', y[i]
    
    print
    h = symbols('g_1:%d'%d)
    g,theta = DSphere(d).gThetaTheta()
    print 'Metric of unit sphere: (gij = gi delta_ij)'
    for i in range(d-1):
        print h[i], ' = ', g[i]
    print
    print 'Area of sphere of radius r:'
    print 'A(S(r)) = ', DSphere(d).AreaDm1(Symbol('r'))
    print
    print 'Volume of compact sphere (ball) of radius r:'
    print 'V(B(r)) = ', DSphere(d).VolD(Symbol('r'))

# =============================================
