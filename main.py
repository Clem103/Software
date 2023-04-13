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
import matplotlib.patches as patches
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)

## Environnement parameters ##
g0 = 9.80665    # Constant for now, might be computed with more precision later
earth_radius = 6.378e6   # In meters, initially the rocket is on the ground



## Parameters of the rocket (Ariane 5G data from Ariane 5 User Manual and Nasa technical Data Sheet) ##
total_mass = np.array([268e3*2,170e3, 10.9e3])   # Initial total mass  of each stage in kg
prop_mass = np.array([237e3*2, 158e3, 9.7e3])
struct_mass = np.copy(total_mass-prop_mass)    # Structure mass in kg
thrust = 9.81*np.array([509.9e3*2, 112.167e3, 2.957e3])     # Thrust in N
Isp = np.array([274.5, 431.2, 321])   # Isp in seconds
V_eq = g0*Isp
m_dot = np.copy(thrust/V_eq)   # Mass flow rate in kg/s
# n_boosters = 2  # Number of booster considered
max_n_stage = 2   # Number of stages considered (0 to max_n_stage)
burn_time = np.sum(prop_mass/m_dot)   # Burn time in seconds


## Simulation parameters ##
end_time = 50000     # End time in seconds
dt = 0.1          # Time step for integrations methods



def rk7():
    pass

# Cf Martin J. L. Turner Rocket Spacecraft prop elements for equations
# Need to find equations for each axis by hand

def compute_accel(r, theta, g0, current_stage, current_step):
    """
    Computes the acceleration of the rocket (assumes earth is flat : the x and y axis are considered fixed)"""

    # Numpy rotation matrix from mobile frame to fixed frame with accel_fixed_frame = rotation_matrix @ accel_mobile_frame
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                                    [np.sin(theta), np.cos(theta)]])

    if prop_mass[current_stage]>=0:
        # In the mobile frame of reference of the rocket

        # Acceleration due to thrust and gravity only (no air drag)
        accel_r = thrust[current_stage]/np.sum(total_mass) * np.sin(pitch_angle_tab[current_step]) - g0
        accel_theta = thrust[current_stage]/np.sum(total_mass) * np.cos(pitch_angle_tab[current_step])

        accel_mobile_frame = np.array([accel_r, accel_theta])
        accel_fixed_frame = rotation_matrix @ accel_mobile_frame
    
    else:

        accel_r = -g0
        accel_theta = 0

        accel_mobile_frame = np.array([accel_r, accel_theta])
        accel_fixed_frame = rotation_matrix @ accel_mobile_frame

    return accel_fixed_frame


def compute_r(x,y):
    return np.sqrt(x**2+y**2)

def compute_alt(r):
    return r - earth_radius

def compute_theta(x,y):
    return np.arctan2(y,x)

def compute_speed(v, a, dt):
    return v + a*dt

def compute_pos(pos, v, dt):
    return pos + dt*v

def stage_sep():
    print("stage sep")
    struct_mass[n_stage] = 0  


time_tab = np.linspace(0,end_time, int(end_time/dt))
pos_tab = []
v_tab = []
M_tab = []
alt_tab = []
v=np.array([0,0])
pos=np.array([0,earth_radius])

pitch_angle_tab1 = np.linspace(np.pi/2,np.pi,int(np.floor(32.43*burn_time*dt)))
pitch_angle_tab = np.append(pitch_angle_tab1, np.pi*(np.ones(int(len(time_tab)-len(pitch_angle_tab1)))))
n_stage = 0


for k in range(len(time_tab)):
    if time_tab[k] == 0:
        pos_tab.append(pos)
        v_tab.append(np.sqrt(v[0]**2+v[1]**2))
        M_tab.append(np.sum(total_mass))
        alt_tab.append(compute_alt(compute_r(pos[0], pos[1])))

    else:
        r = compute_r(pos[0], pos[1])
        
        if r < earth_radius:
            print("Crash")
            break
        
        alt= compute_alt(r)
        theta = compute_theta(pos[0], pos[1])
        a = compute_accel(r, theta, g0, n_stage, k)
        v = compute_speed(v, a, dt)
        pos = compute_pos(pos, v, dt)

        v_tab.append(np.sqrt(v[0]**2+v[1]**2))
        pos_tab.append(pos)
        alt_tab.append(alt)
        prop_mass[n_stage] -= m_dot[n_stage]*dt

        if prop_mass[n_stage]<= 0:
            if n_stage<max_n_stage:
                stage_sep()
                n_stage += 1
            else:
                prop_mass[n_stage] = 0
        
        total_mass = prop_mass + struct_mass
        M_tab.append(np.sum(total_mass))

        

pos_tab = np.array(pos_tab)
v_tab = np.array(v_tab)

fig, axes = plt.subplots(2,2)
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

ax1.plot(pos_tab[:,0], pos_tab[:,1], label='Position', color='r')
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.set_title('Rocket trajectory')
ax1.grid()
ax1.legend()
circle1 = plt.Circle((0, 0), earth_radius, color='b', alpha=0.5)
ax1.axis('equal')
ax1.add_patch(circle1)
ax2.plot(time_tab[:len(v_tab)], v_tab, label='Speed')
ax2.legend()
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Speed (m/s)')
ax2.set_title('Rocket speed')
ax2.grid()
ax3.plot(time_tab[:len(v_tab)], alt_tab, label='Altitude')
ax3.legend()
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Altitude (m)')
ax3.set_title('Rocket altitude')
ax3.grid()
ax4.plot(time_tab[:len(v_tab)], pitch_angle_tab[:len(v_tab)], label='Pitch angle')
ax4.legend()
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('Pitch angle (rad)')
ax4.set_title('Pitch angle')
ax4.grid()
plt.show()
