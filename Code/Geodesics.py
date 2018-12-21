from __init__ import *

# Extra forces set to zero for this run
# -------------------------------------
def XFORCE(t, x, v):
    return np.zeros(len(x)+1, float)

# Input metric and coordinate system
# ----------------------------------
D = 4
d = D-1
#M,J,Q = 1.0, 0.0, 0.0
M,J,Q = 1.0, 0.9, 0.0
#M,J,Q = 1.0, 0.9, 0.3

g,x = metric(M,J,Q)  #g,x = metric(D,M,Q)

# Initial Conditions
# ------------------
rSigma = 1.e-2
t0 = 0.0
x0 = np.zeros(d, float)
v0 = np.zeros(d, float)

#x0[0] = 2.0/rSigma
x0[0] = 1.72223189
x0[d-1] = 1.15216234
for i in range(1,d-1):
    x0[i] = np.pi/2
v0[d-1] = 0.09949532
#v0[0] = -1.e-8

q0 = g_ginv_Gamma( x, g )
q = GEODESICS_SOLVER( q0.gf, q0.dgf, q0.Gammaf, XFORCE )
q.SetInitialConditions( t0, x0, v0 )
q.SetXForceArgs( *[] )


# FALLING PARTICLE FOR D=4 SCHWARZSCHILD
#x0[0] = 3.0
#for i in range(1,d-1):
#    x0[i] = np.pi/2
#v0[0] = -0.3
#v0[d-1] = 0.03

# REISNET-NORDSTROM INITIAL CONDITIONS FOR BOUNDED TRAJECTORY WITH T~240
#x0[0] = 10.0
#for i in range(1,d-1):
#    x0[i] = np.pi/2
#v0[0] = -0.03
#v0[d-1] = 0.03

print 'Initial Conditions'
print '------------------'
print 't_0 = ', q.t0
print 'x_0 = ', q.x0
print 'v_0 = ', q.v0
tN = float(input('Final time instant: tN = '))

INTEGRATOR, INTEGRATOR_LABEL, INTEGRATOR_ERROR_ORDER, SYMPLECTIC_FLAG, ADAPTIVE_FLAG = getIntegratorRUN()
q.SetIntegrator( INTEGRATOR, INTEGRATOR_ERROR_ORDER )

# Solve geodesics
# ---------------
if ADAPTIVE_FLAG == True:
    sigma = float(input('Tolerance per step: sigma = '))
    print sigma
    dt0 = sigma*(tN-t0)
    q.evolveAdaptive( tN, sigma, dt0 )
else:
    N = int(input('Number of time steps: N = '))
    print N
    q.evolve( tN, N, SYMPLECTIC_FLAG )

NH = 10000
rHtmp = eventHorizon(M,J,Q) #rHtmp = eventHorizon(D,M,J,Q)
rH = rHtmp*np.ones(NH, float)
thetaH = np.pi/2*np.ones(NH, float)  #thetaH = np.pi*np.arange(NH)/NH
phiH = 2*np.pi*np.arange(NH)/NH
xH = rH*np.sin(thetaH)*np.cos(phiH)
yH = rH*np.sin(thetaH)*np.sin(phiH)
zH = rH*np.cos(thetaH)

fig, ax = plt.subplots(d,2)
for i in range(d):
    ax[i,0].plot(q.t, q.x[i,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    if i==0:
        ax[i,0].plot(q.t, rH[:q.N],ls='--',color='r',label='Event Horizon')
    ax[i,0].set_ylabel(r'${}(t)$'.format(x[i+1]))
    ax[i,0].legend(loc='upper right')
    ax[i,1].plot(q.t, q.v[i,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax[i,1].set_ylabel(r'$v_{%r}(t)$'%x[i+1])
    ax[i,1].yaxis.tick_right()
    ax[i,1].yaxis.set_label_position('right')
    ax[i,1].legend(loc='upper right')
ax[d-1,0].set_xlabel(r'$t$')
ax[d-1,1].set_xlabel(r'$t$')
plt.savefig("xv_t.pdf",bbox_inches="tight")

U = q.FourMomenta()

fig3, ax3 = plt.subplots(D,1)
for mu in range(D):
    ax3[mu].plot(q.t, U[mu,:],ls='',marker='.',ms=3.0,label=INTEGRATOR_LABEL)
    ax3[mu].set_ylabel(r'$p_%d(t)/m$'%mu)
    ax3[mu].legend(loc='upper right')
ax3[D-1].set_xlabel(r'$t$')
plt.savefig("U_t.pdf",bbox_inches="tight")

t = q.t
x = q.x[0,:]*np.sin(q.x[1,:])*np.cos(q.x[2,:])
y = q.x[0,:]*np.sin(q.x[1,:])*np.sin(q.x[2,:])
z = q.x[0,:]*np.cos(q.x[1,:])

# Animation production... Takes too long!
# ---------------------------------------
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#fig1, ax1 = plt.subplots()

#def animate(n):
##    ax1.plot(x[n],y[n],marker='o',ms=1.0,color='g')
#    ax1.plot(x[n],y[n],color=(colors[n],0.1,colors[frms-1-n]),marker='o',ms=1.0)
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

datafile = open('t_x.txt','w+')
for n in range(q.N):
    datafile.write('%f    %f   %f   %f\n'%(q.t[n],q.x[0,n],q.x[1,n],q.x[2,n]))

datafile = open('t_v.txt','w+')
for n in range(q.N):
    datafile.write('%f    %f   %f   %f\n'%(q.t[n],q.v[0,n],q.v[1,n],q.v[2,n]))

datafile = open('t_U.txt','w+')
for n in range(q.N):
    datafile.write('%f    %f   %f   %f   %f\n'%(q.t[n],U[0,n],U[1,n],U[2,n],U[3,n]))

plt.show()
