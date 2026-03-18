from aero_main import TakeoffCoeff, CruiseCoeff
import matplotlib.pyplot as plt
import numpy as np

# drag polar??
C_D = CruiseCoeff.CD_tot
C_L = CruiseCoeff.CL

plt.figure()
plt.plot(C_D, C_L, label = 'cruise')
plt.xlabel('$C_D$')
plt.ylabel('$C_L$')
plt.legend()
plt.show()


plt.figure()
plt.plot(C_D, C_L, label = 'takeoff')
plt.xlabel('$C_D$')
plt.ylabel('$C_L$')
plt.legend()
plt.show()