# Programmed using Python 3.10
# Author : Clément Vellu

# Package requirements : 
# - Numpy 1.24.2
# - Matplotlib 3.7.1

import numpy as np
import matplotlib.pyplot as plt

## Environnement parameters ##
g0 = 9.80665    # Constant for now, might be computed with more precision later



## Parameters of the rocket (Ariane 5 data) ##
M = 746000   # Initial mass in kg
Isp = 274   # Isp in seconds
m_dot = 2000   # Mass flow rate in kg/s
n_boosters = 2  # Number of booster considered


## Simulation parameters ##
end_time = 1000     # End time in seconds
dt = 0.01           # Time step for integrations methods



def rk7():
    pass

def compute_accel(m_dot, Isp, g0, M):
    if M > 0:
        return g0*(n_boosters*Isp*m_dot/M - 1)
    else:
        return 0

def compute_speed(v, a, dt):
    return v + a*dt

def compute_pos(x, v, dt):
    return x + dt*v

time_tab = np.linspace(0,end_time, int(end_time/dt))
pos_tab = []
v_tab = []
M_tab = []
v=0
x=0


for t in time_tab:
    if t == 0:
        pos_tab.append(x)
        v_tab.append(v)
        M_tab.append(M)

    else:
        a = compute_accel(m_dot, Isp, g0, M)
        v = compute_speed(v, a, dt)
        x = compute_pos(x, v, dt)

        v_tab.append(v)
        pos_tab.append(x)
        M -= m_dot*dt
        M_tab.append(M)

plt.plot(time_tab, pos_tab, label='Position')
plt.legend()
plt.figure()
plt.plot(time_tab, v_tab, label='Speed')
plt.legend()
plt.figure()
plt.plot(time_tab, M_tab, label = 'Mass')
plt.legend()
plt.show()


