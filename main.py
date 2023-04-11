# Programmed using Python 3.10
# Author : Clément Vellu

# Package requirements : 
# - Numpy 1.24.2
# - Matplotlib 3.7.1

# What needs to be done : 
# - Implement trajectory equation with linear pitch angle
# - Implement air drag
# - Implement other rocket type

import numpy as np
import matplotlib.pyplot as plt

## Environnement parameters ##
g0 = 9.80665    # Constant for now, might be computed with more precision later



## Parameters of the rocket (Ariane 5G data from Ariane 5 User Manual and Nasa technical Data Sheet) ##
total_mass = np.array([268e3*2,170e3, 10.9e3])   # Initial total mass  of each stage in kg
prop_mass = np.array([237e3*2, 158e3, 9.7e3])
struct_mass = np.copy(total_mass-prop_mass)    # Structure mass in kg
thrust = 9.81*np.array([509.9e3*2, 112.167e3, 2.957e3])     # Thrust in N
Isp = np.array([274.5, 431.2, 321])   # Isp in seconds
V_eq = g0*Isp
m_dot = np.copy(thrust/V_eq)   # Mass flow rate in kg/s
# n_boosters = 2  # Number of booster considered
max_n_stage = 2


## Simulation parameters ##
end_time = 1800     # End time in seconds
dt = 0.01          # Time step for integrations methods



def rk7():
    pass

# Cf Martin J. L. Turner Rocket Spacecraft prop elements for equations
# Need to find equations for each axis by hand

def compute_accel(m_dot, Isp, g0, M, current_stage, current_step):
    """
    Computes the acceleration of the rocket (assumes earth is flat : the x and y axis are considered fixed)"""
    accel_vect_x = 0
    accel_vect_y = 0
    if prop_mass[current_stage]>=0:
        # return g0*(n_boosters*Isp*m_dot/M - 1) Legacy formula
        accel_vect_x = thrust[current_stage]/np.sum(total_mass) * np.cos(pitch_angle_tab[current_step])
        accel_vect_y = thrust[current_stage]/np.sum(total_mass) * np.sin(pitch_angle_tab[current_step]) - g0
    else:
        accel_vect_y = -g0

    return np.array([accel_vect_x, accel_vect_y])

def compute_speed(v, a, dt):
    return v + a*dt

def compute_pos(x, v, dt):
    return x + dt*v

def stage_sep():
    print("stage sep")
    struct_mass[n_stage] = 0

time_tab = np.linspace(0,end_time, int(end_time/dt))
pos_tab = []
v_tab = []
M_tab = []
v=np.array([0,0])
pos=np.array([0,0])

pitch_angle_tab = np.linspace(np.pi/2,0,len(time_tab))

# struct_mass = struct_mass_s1
# n_boosters = n_boosters_s1
# m_dot = m_dot_s1
# Isp = Isp_s1
n_stage = 0


for k in range(len(time_tab)):
    if time_tab[k] == 0:
        pos_tab.append(pos)
        v_tab.append(v)
        M_tab.append(np.sum(total_mass))

    else:
        a = compute_accel(m_dot, Isp, g0, total_mass, n_stage, k)
        v = compute_speed(v, a, dt)
        pos = compute_pos(pos, v, dt)

        v_tab.append(v)
        pos_tab.append(pos)
        prop_mass[n_stage] -= m_dot[n_stage]*dt
        total_mass = prop_mass + struct_mass
        M_tab.append(np.sum(total_mass))

        if prop_mass[n_stage]<= 0:
            if n_stage<max_n_stage:
                stage_sep()
                n_stage += 1

pos_tab = np.array(pos_tab)
v_tab = np.array(v_tab)

plt.plot(pos_tab[:,0], pos_tab[:,1], label='Position')
plt.legend()
plt.figure()
plt.plot(v_tab[:,0], v_tab[:,1], label='Speed')
plt.legend()
plt.figure()
plt.plot(time_tab, M_tab, label = 'Mass')
plt.legend()
plt.show()


