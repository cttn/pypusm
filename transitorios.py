import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from SisEqDif import pusm_re as pusm

t_ini = 0   #dpp
t_fin = 100  #dpp
numpoints = 1000

t = [t_ini + (t_fin - t_ini) * float(i) / (numpoints - 1) for i in range(numpoints)]

concs_ini = [1.0, 1.0, 1.0, 1.0]

def FPP(tiempo):
    if (tiempo > 0.0 and tiempo < 20.0):
        return 1.0
    elif (tiempo > 20.0 and tiempo < 40.0):
        return 0.5
    elif (tiempo > 40.0 and tiempo < 60.0):
        return 1.0
    else:
        return 0.0

soluciones = odeint(pusm, concs_ini, t, args=(FPP,), hmax=0.01)

Esm, Epm, Enp, Epu = soluciones.T #Concentraciones referidas a 1FPP
potencia = [FPP(i) for i in t]    #Vector de potencias

#Cambios en la reactividad para combustible en equilibrio
coef_Pu = 270.6 #mk
coef_Sm =  -4.6 #mk

Drho = np.zeros(numpoints)
for i in range(len(t)):
    Drho[i] = (Esm[i]-1.0)*coef_Sm + (Epu[i]-1.0)*coef_Pu
