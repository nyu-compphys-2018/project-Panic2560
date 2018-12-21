import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from METRIC.metric import *
from SOLVER.solver import *
from SOLVER.integrators import *
from SOLVER.g_ginv_Gamma import *
from SOLVER.PenroseSplitIC import *

def getIntegratorRUN():
    INTEGRATORS_LABELS = ( 'Forward Euler', 'RK2', 'RK4', 'Leapfrog', 'Verlet', 'Adaptive Forward Euler', 'Adaptive RK2', 'Adaptive RK4')
    print "Choose Integrator:"
    for i in range(8):
        print '[',i,'] ', INTEGRATORS_LABELS[i]
    iINT = int(input('Integrator ID = '))
    print iINT
    if iINT == 0:
        INTEGRATOR = F_EULER
        INTEGRATOR_ERROR_ORDER = 1
        SYMPLECTIC_FLAG = False
        ADAPTIVE_FLAG = False
    elif iINT == 1:
        INTEGRATOR = RK2
        INTEGRATOR_ERROR_ORDER = 3
        SYMPLECTIC_FLAG = False
        ADAPTIVE_FLAG = False
    elif iINT == 2:
        INTEGRATOR = RK4
        INTEGRATOR_ERROR_ORDER = 5
        SYMPLECTIC_FLAG = False
        ADAPTIVE_FLAG = False
    elif iINT == 3:
        INTEGRATOR = LEAPFROG
        INTEGRATOR_ERROR_ORDER = 2
        SYMPLECTIC_FLAG = True
        ADAPTIVE_FLAG = False
    elif iINT == 4:
        INTEGRATOR = VERLET
        INTEGRATOR_ERROR_ORDER = 2
        SYMPLECTIC_FLAG = True
        ADAPTIVE_FLAG = False
    elif iINT == 5:
        print 'FOOL YA! Forward Euler cannot be implemented with an adaptive step because it is a first order integrator.'
        print
        getIntegratorRUN()
    elif iINT == 6:
        INTEGRATOR = RK2
        INTEGRATOR_ERROR_ORDER = 3
        SYMPLECTIC_FLAG = False
        ADAPTIVE_FLAG = True
    elif iINT == 7:
        INTEGRATOR = RK4
        INTEGRATOR_ERROR_ORDER = 5
        SYMPLECTIC_FLAG = False
        ADAPTIVE_FLAG = True
    else:
        print 'INVALID VALUE! CHOOSE BETWEEN 0 - 7'
        print
        getIntegratorRUN()

    return INTEGRATOR, INTEGRATORS_LABELS[iINT], INTEGRATOR_ERROR_ORDER, SYMPLECTIC_FLAG, ADAPTIVE_FLAG
