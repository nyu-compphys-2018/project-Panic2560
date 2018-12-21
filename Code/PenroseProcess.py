from __init__ import *

# Extra forces set to zero for this run
# -------------------------------------
def XFORCE(t, x, v):
    return np.zeros(len(x)+1, float)

# Input metric and coordinate system
# ----------------------------------
D = 4
d = D-1
#M, J, Q = 1.0, 0.9, 0.2
M, J, Q = 1.0, 0.9, 0.0

g,x = metric(M,J,Q)

# Initial Conditions for particle 1
# ---------------------------------
t01 = 3176.53500
x01 = np.zeros(d, float)
v01 = np.zeros(d, float)
x01[0] = 6.103719
for i in range(1,d-1):
    x01[i] = np.pi/2
x01[d-1] = 0.072048
v01[0] = -0.386224
v01[d-1] = 0.007694

q0 = g_ginv_Gamma( x, g )
q1 = GEODESICS_SOLVER( q0.gf, q0.dgf, q0.Gammaf, XFORCE )
q1.SetInitialConditions( t01, x01, v01 )
q1.v0i_umuConstaint( d-1, d, 0.0 )
q1.SetXForceArgs( *[] )

print 'Initial Conditions for first process'
print '------------------------------------'
print 't_01 = ', q1.t0
print 'x_01 = ', q1.x0
print 'v_01 = ', q1.v0
#tN = float(input('Final time instant for first process: tN1 = '))
#print tN1
tN1 = 3191.457500 # r = 1.720577, phi = 1.156470, vr = -0.084983, vd = 0.222024
print 'Final time instant for first process: tN1 = ', tN1

INTEGRATOR, INTEGRATOR_LABEL, INTEGRATOR_ERROR_ORDER, SYMPLECTIC_FLAG, ADAPTIVE_FLAG = getIntegratorRUN()
q1.SetIntegrator( INTEGRATOR, INTEGRATOR_ERROR_ORDER )

# Solve geodesics for particle 1
# ------------------------------
if ADAPTIVE_FLAG == True:
    sigma1 = float(input('Tolerance per step for first process: sigma1 = '))
    print sigma1
    dt01 = sigma1*(tN1-t01)
    q1.evolveAdaptive( tN1, sigma1, dt01 )
else:
#    N1 = int(input('Number of time steps for first process: N1 = '))
#    print N1
    N1 = 1000
    print 'Number of time steps for first process: N1 = ', N1
    q1.evolve( tN1, N1, SYMPLECTIC_FLAG )

print
print '==========================================================='
print '==========================================================='
print

# Initial positions for particles 2 and 3
# ---------------------------------------
x02 = np.zeros(d, float)
x03 = np.zeros(d, float)
t02 = t03 = q1.t[q1.N-1]
for i in range(d):
    x02[i] = q1.x[i,q1.N-1]
    x03[i] = q1.x[i,q1.N-1]


y1 = np.zeros(D, float)
y1[0] = q1.t[q1.N-1]
y1[1:] = q1.x[:,q1.N-1]
g0 = q0.gf(*y1)

m1 = 1.0

v02 = np.zeros(d, float)
v03 = np.zeros(d, float)
v03[0] = 0.0
v03[1] = 0.0
omega = -g0[0][0]/g0[0][d]
#psi = np.sqrt(-g0[0][0]*g0[d][d]+(g0[0][d])**2)/g0[d][d]
#print -psi, +psi
#v3[2] = omega - 1.e-6                                                      
v03[2] = omega - 5.4820298e-2   # Set v03[2] such that E3 < 0

q = PENROSE_SPLIT_IC( g0, m1, q1.v[:,q1.N-1] )
q.Setv3( v03 )
q.m2m3v2d( )    # Solve to find m2/m1, m3/m1 and v20
v02 = q.v2

# Solve geodesics for particles 2 and 3
# -------------------------------------
q2 = GEODESICS_SOLVER( q0.gf, q0.dgf, q0.Gammaf, XFORCE )
q2.SetInitialConditions( t02, x02, v02 )
q2.SetXForceArgs( *[] )
q2.SetIntegrator( INTEGRATOR, INTEGRATOR_ERROR_ORDER )

q3 = GEODESICS_SOLVER( q0.gf, q0.dgf, q0.Gammaf, XFORCE )
q3.SetInitialConditions( t03, x03, v03 )
q3.SetXForceArgs( *[] )
q3.SetIntegrator( INTEGRATOR, INTEGRATOR_ERROR_ORDER )

print 'Initial Conditions for second process (Outgoing particle)'
print '---------------------------------------------------------'
print 't_02 = ', q2.t0
print 'x_02 = ', q2.x0
print 'v_02 = ', q2.v0
#tN2 = float(input('Final time instant for second process (Outgoing particle): tN2 = '))
#print tN2
tN2 = 3500.0
print 'Final time instant for second process (Outgoing particle): tN2 = ', tN2
print
print 'Initial Conditions for second process (Ingoing particle)'
print '--------------------------------------------------------'
print 't_03 = ', q3.t0
print 'x_03 = ', q3.x0
print 'v_03 = ', q3.v0
#tN3 = float(input('Final time instant for second process (Ingoing particle): tN3 = '))
#print tN3
tN3 = q1.t[q1.N-1] + 15.0
print 'Final time instant for second process (Ingoing particle): tN3 = ', tN3

if ADAPTIVE_FLAG == True:
    sigma2 = float(input('Tolerance per step for second process (Outgoing particle): sigma2 = '))
    print sigma2
    sigma3 = float(input('Tolerance per step for second process (Ingoing particle): sigma3 = '))
    print sigma3
    dt02 = sigma2*(tN2-t02)
    dt03 = sigma3*(tN3-t03)
    q2.evolveAdaptive( tN2, sigma2, dt02 )
    q3.evolveAdaptive( tN3, sigma3, dt03 )
else:
#    N2 = int(input('Number of time steps for second process: N2 = '))
#    print N2
#    N3 = int(input('Number of time steps for second process: N3 = '))
#    print N3
    N2 = 20000
    print 'Number of time steps for first process: N2 = ', N2
    N3 = 20000
    print 'Number of time steps for first process: N3 = ', N3
    q2.evolve( tN2, N2, SYMPLECTIC_FLAG )
    q3.evolve( tN3, N3, SYMPLECTIC_FLAG )

N1, N2, N3 = q1.N, q2.N, q3.N
N23 = max(N2,N3)
t = np.zeros(N1+N23, float)
rH = eventHorizon(M,J,Q)*np.ones(N1+N23,float)
rERGO = np.zeros(N1+N23, float)
for n in range(N1):
    t[n] = q1.t[n]
    rERGO[n] = ergoSurface(q1.x[1,n],M,J,Q)
for n in range(N23):
    if tN2>tN3:
        t[N1+n] = q2.t[n]
        rERGO[N1+n] = ergoSurface(q2.x[1,n],M,J,Q)
    else:
        t[N1+n] = q3.t[n]
        rERGO[N1+n] = ergoSurface(q3.x[1,n],M,J,Q)

fig, ax = plt.subplots(d,2)
for i in range(d):
    ax[i,0].plot(q1.t, q1.x[i,:],ls='',marker='.',ms=3.0)
    ax[i,0].plot(q2.t[:800], q2.x[i,:800],ls='',marker='.',ms=3.0)
    ax[i,0].plot(q3.t, q3.x[i,:],ls='',marker='.',ms=3.0)
    if i==0:
        ax[i,0].plot(t,rH,ls='--',color='r',label='Event Horizon')
        ax[i,0].plot(t,rERGO,ls='--',color='g',label='Ergosurface')
    ax[i,0].set_ylabel(r'${}(t)$'.format(x[i+1]))
    ax[i,0].legend(loc='upper right')
    ax[i,1].plot(q1.t, q1.v[i,:],ls='',marker='.',ms=3.0)
    ax[i,1].plot(q2.t[:800], q2.v[i,:800],ls='',marker='.',ms=3.0)
    ax[i,1].plot(q3.t, q3.v[i,:],ls='',marker='.',ms=3.0)
    ax[i,1].set_ylabel(r'$v_{%r}(t)$'%x[i+1])
    ax[i,1].yaxis.tick_right()
    ax[i,1].yaxis.set_label_position('right')
ax[d-1,0].set_xlabel(r'$t$')
ax[d-1,1].set_xlabel(r'$t$')
plt.savefig("t_x_v.pdf",bbox_inches="tight")

U1 = q.m1*q1.FourMomenta()
U2 = q.m2*q2.FourMomenta()
U3 = q.m3*q3.FourMomenta()

fig3, ax3 = plt.subplots(D,1)
for mu in range(D):
    ax3[mu].plot(q1.t, U1[mu,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax3[mu].plot(q2.t, U2[mu,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax3[mu].plot(q3.t[:N2], U3[mu,:N2],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax3[mu].plot(q2.t[:N2], U3[mu,:N2]+U2[mu,:],ls='--')
    ax3[mu].set_ylabel(r'$p_%d(t)$'%mu)
    ax3[mu].legend(loc='upper right')
ax3[D-1].set_xlabel(r'$t$')
plt.savefig("U_t.pdf",bbox_inches="tight")

fig4, ax4 = plt.subplots(D,1)
for mu in range(D):
    ax4[mu].plot(q1.t, U1[mu,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax4[mu].plot(q2.t, U2[mu,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax4[mu].plot(q3.t, U3[mu,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax4[mu].set_ylabel(r'$p_%d(t)/m$'%mu)
    ax4[mu].legend(loc='upper right')
ax4[D-1].set_xlabel(r'$t$')
plt.savefig("U_t_Full.pdf",bbox_inches="tight")

# Animation production... Takes too long!
# ---------------------------------------
#t = q.t
#x = q.x[0,:]*np.sin(q.x[1,:])*np.cos(q.x[2,:])
#y = q.x[0,:]*np.sin(q.x[1,:])*np.sin(q.x[2,:])
#z = q.x[0,:]*np.cos(q.x[1,:])

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#fig1, ax1 = plt.subplots()

#def animate(n):
#    ax1.plot(x[n],y[n],marker='o',ms=1.0,color='g')
##    ax1.plot(x[n],y[n],color=(colors[n],0.1,colors[frms-1-n]),marker='o',ms=1.0)
#    plt.legend(['t = %f'%t[n]],loc=9)

#frms = q.N

#colors = np.arange(frms)/float(frms)
#lim = max(np.max(x),np.max(y))
#ax1.set_facecolor('black')
#ax1.set_xlim(-1.1*lim,1.1*lim)
#ax1.set_ylim(-1.1*lim,1.1*lim)
#ax1.plot(xH,yH,'r--')

#anim = animation.FuncAnimation(fig1, animate, frames=frms-1, interval=10)
#anim.save('Rotating.mp4',writer=writer)

datafile1x = open('P1_x.txt','w+')
datafile1v = open('P1_v.txt','w+')
datafile1U = open('P1_U.txt','w+')
for n in range(N1):
    datafile1x.write('%f   %f   %f   %f\n'%(q1.t[n],q1.x[0,n],q1.x[1,n],q1.x[2,n]))
    datafile1v.write('%f   %f   %f   %f\n'%(q1.t[n],q1.v[0,n],q1.v[1,n],q1.v[2,n]))
    datafile1U.write('%f   %f   %f   %f   %f\n'%(q1.t[n],U1[0,n],U1[1,n],U1[2,n],U1[3,n]))

datafile2x = open('P2_x.txt','w+')
datafile2v = open('P2_v.txt','w+')
datafile2U = open('P2_U.txt','w+')
for n in range(N2):
    datafile2x.write('%f   %f   %f   %f\n'%(q2.t[n],q2.x[0,n],q2.x[1,n],q2.x[2,n]))
    datafile2v.write('%f   %f   %f   %f\n'%(q2.t[n],q2.v[0,n],q2.v[1,n],q2.v[2,n]))
    datafile2U.write('%f   %f   %f   %f   %f\n'%(q2.t[n],U2[0,n],U2[1,n],U2[2,n],U2[3,n]))

datafile3x = open('P3_x.txt','w+')
datafile3v = open('P3_v.txt','w+')
datafile3U = open('P3_U.txt','w+')
for n in range(N3):
    datafile3x.write('%f   %f   %f   %f\n'%(q3.t[n],q3.x[0,n],q3.x[1,n],q3.x[2,n]))
    datafile3v.write('%f   %f   %f   %f\n'%(q3.t[n],q3.v[0,n],q3.v[1,n],q3.v[2,n]))
    datafile3U.write('%f   %f   %f   %f   %f\n'%(q3.t[n],U3[0,n],U3[1,n],U3[2,n],U3[3,n]))

plt.show()
