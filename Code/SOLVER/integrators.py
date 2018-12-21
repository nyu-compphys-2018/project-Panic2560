from __init__ import *
# ====================================================================

# ====================================================================
# =================== NON-SYMPLECTIC INTEGRATORS =====================
# ====================================================================
def F_EULER( t, x, v, dt, FORCE, *Fargs ):

    xn = x + dt*v
    vn = v + dt*FORCE( t, x, v, *Fargs )

    return xn,vn
# ====================================================================
# ====================================================================
def RK2( t, x, v, dt, FORCE, *Fargs ):

    k1x = dt*v
    k1v = dt*FORCE( t, x, v, *Fargs )
    k2x = dt*(v + k1v/2)
    k2v = dt*FORCE( t+dt/2, x+k1x/2, v+k1v/2, *Fargs )

    xn = x + k2x
    vn = v + k2v

    return xn,vn
# ====================================================================
# ====================================================================
def RK4( t, x, v, dt, FORCE, *Fargs ):

    k1x = dt*v
    k1v = dt*FORCE( t, x, v, *Fargs )
    k2x = dt*(v + k1v/2)
    k2v = dt*FORCE( t+dt/2, x+k1x/2, v+k1v/2, *Fargs )
    k3x = dt*(v + k2v/2)
    k3v = dt*FORCE( t+dt/2, x+k2x/2, v+k2v/2, *Fargs )
    k4x = dt*(v + k3v)
    k4v = dt*FORCE( t+dt, x+k3x, v+k3v, *Fargs )

    xn = x + (k1x + 2*k2x + 2*k3x + k4x)/6
    vn = v + (k1v + 2*k2v + 2*k3v + k4v)/6

    return xn,vn
# ====================================================================

# ====================================================================
# ===================== SYMPLECTIC INTEGRATORS =======================
# ====================================================================
def LEAPFROG( t0, x0, v0, tN, N, FORCE, *Fargs ):
    d = len(x0)

    dt = (tN - t0)/float(N)
    t = t0 + dt*np.arange(N)
    x = np.zeros([d,N], float)
    v = np.zeros([d,N], float)

    t[0] = t0
    x[:,0] = x0
    v[:,0] = v0

    xaux = x0 + (dt/2)*v0
    vaux = v0 + (dt/2)*FORCE( t0, x0, v0, *Fargs )

    for n in range(N-1):
        x[:,n+1] = x[:,n] + dt*vaux
        v[:,n+1] = v[:,n] + dt*FORCE( t[n]+dt/2, xaux, vaux, *Fargs )
        xaux += dt*v[:,n+1]
        vaux += dt*FORCE( t[n+1], x[:,n+1], v[:,n+1], *Fargs )
    # ----- GRAGG'S STEP ------
    x[:,N-1] = 0.5*(x[:,N-1] + x[:,N-2] + dt*v[:,N-1])
    v[:,N-1] = 0.5*(v[:,N-1] + v[:,N-2] + dt*FORCE( t[N-1], x[:,N-1], v[:,N-1], *Fargs ))
    # -------------------------

    return t, x, v
# ====================================================================
# ====================================================================
def VERLET( t0, x0, v0, tN, N, FORCE, *Fargs ):
    # Here, the FORCE does not depend on the velocity \vec{v} but it is
    # kept explicitly in the arguments of the FORCE for it to be callable
    # in the SOLVE_NEWTON class
    d = len(x0)

    dt = (tN - t0)/float(N)
    t = t0 + dt*np.arange(N)
    x = np.zeros([d,N], float)
    v = np.zeros([d,N], float)

    t[0] = t0
    x[:,0] = x0
    v[:,0] = v0

    vaux = v0 + (dt/2)*FORCE( t0, x0, v0, *Fargs )

    for n in range(N-1):
        x[:,n+1] = x[:,n] + dt*vaux
        kv = dt*FORCE( t[n+1], x[:,n+1], v0, *Fargs )
        v[:,n+1] = v[:,n] + kv/2
        vaux += kv
    # ----- GRAGG'S STEP ------
    x[:,N-1] = 0.5*(x[:,N-1] + x[:,N-2] + dt*v[:,N-1])
    v[:,N-1] = 0.5*(v[:,N-1] + v[:,N-2] + dt*FORCE( t[N-1], x[:,N-1], v[:,N-1], *Fargs ))
    # -------------------------

    return t, x, v
