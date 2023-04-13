# Name : Rocket trajectory simulator
# Current version number : 4.0
# Author : Clément Vellu
# Initial coding : 20/03/2023
# Description : This program simulates the trajectory of a rocket lauching from surface of Earth using Euler's integration scheme. 
#               The rocket modelled is Ariane 5 rocket with 3 stages
#               The simulator can be used in three modes : linearly increasing pitch angle, gravity turn, constant pitch angle
# Algorithm used : Euler integration scheme (coded from scratch)
# Expected working range of parameters :
#      - Rocket parameters : Positive values for all parameters
#      - Simulation parameters :
#           - End time : Positive value
#           - Time step : Positive value
#           - Pitch slope time : Positive value
#           - Initial pitch difference : Any value
#           - Constant pitch angle : Any value
#           - User input : 0, 1 or 2
# Tested working range of parameters :
#      - Rocket parameters : Ariane 5 rocket parameters
#      - Simulation parameters :
#           - End time : 100/500/1000/3000/5000 seconds
#           - Time step : 0.01/0.1 seconds
#           - Pitch slope time : 100/1000/2000 seconds
#           - Initial pitch difference : 1e-8 radians
#           - Constant pitch angle : -np.pi/8 radians
#           - User input : 0, 1 and 2
# List of inputs : Mode (values 0, 1 or 2)
# List of outputs : Graphical plots of trajectory in the fixed reference frame (x,y), speed over time in the fixed reference frame, 
#                   altitude (distance from the surface of the Earth) over time of the rocket, pitch angle over time (for mode 0 and 2 in the mobile frame of reference, for mode 1 in the fixed reference frame)
# Package requirements and versions :
# - Numpy 1.24.2
# - Matplotlib 3.7.1
# 
# Python version used for development : 3.10
# List of changes and date of release : see https://github.com/Clem103/Software/commits/master/ for a detailed list

import numpy as np
import matplotlib.pyplot as plt

############ Environnement parameters ############

g0 = 9.80665             # In m/s^2
earth_radius = 6.378e6   # In meters, initially the rocket is on the earth surface


############ Parameters of the rocket (Ariane 5G data from Ariane 5 User Manual and Nasa technical Data Sheet) ############
total_mass = np.array([268e3*2,170e3, 10.9e3])               # Initial total mass of each stage in kg
prop_mass = np.array([237e3*2, 158e3, 9.7e3])                # Initial propellant mass of each stage in kg
struct_mass = np.copy(total_mass-prop_mass)                  # Structure mass of each stage in kg
thrust = 9.81*np.array([509.9e3*2, 112.167e3, 2.957e3])      # Thrust for each stage in N
Isp = np.array([274.5, 431.2, 321])                          # Isp for each stage in seconds
V_eq = g0*Isp                                                # Exit velocity for each stage in m/s
m_dot = np.copy(thrust/V_eq)                                 # Mass flow rate for each stage in kg/s
max_n_stage = 2                                              # Number of stages considered (0 to max_n_stage)
burn_time = np.sum(prop_mass/m_dot)                          # Burn time in seconds


############ Simulation parameters ############
# Can be safely modified by the user #

end_time = 3000                 # End time in seconds
dt = 0.1                        # Time step for integration method in seconds
pitch_slope_time = 2000         # MODE 0 : Time in seconds during which the pitch angle is linearly increased from vertical to horizontal (in mobile reference frame) 
init_pitch_diff = 1e-8          # MODE 1 : Initial pitch angle difference between the vertical axis and the rocket in the fixed reference frame (in radians) (used to initiate gravity turn) 
constant_pitch = -np.pi/8       # MODE 2 : Constant pitch angle between the vertical axis and the rocket in the mobile frame of reference (in radians)


def compute_accel_mobile(beta, current_stage, current_step):
    """
    Computes the acceleration of the rocket in the mobile frame of reference (r,beta) and returns it in the fixed reference frame (x,y)
    Inputs : 
        - beta : angle between mobile frame of reference origin and fixed reference frame origin (in radians)
        - current_stage : current stage of the rocket (0 to max_n_stage)
        - current_step : current step of the simulation (0 to end_time/dt)
    Outputs :
        - accel_fixed_frame : vector acceleration of the rocket in the fixed reference frame (x,y) (in m/s^2)
    """

    # Numpy rotation matrix from mobile frame to fixed frame with accel_fixed_frame = rotation_matrix @ accel_mobile_frame
    rotation_matrix = np.array([[np.cos(beta), -np.sin(beta)],
                                [np.sin(beta), np.cos(beta)]])

    if prop_mass[current_stage]>=0:
        # In the mobile frame of reference of the rocket

        # Acceleration due to thrust and gravity only (no air drag)
        accel_r = thrust[current_stage]/np.sum(total_mass) * np.cos(pitch_angle_tab[current_step]) - g0
        accel_beta = thrust[current_stage]/np.sum(total_mass) * np.sin(pitch_angle_tab[current_step])

        accel_mobile_frame = np.array([accel_r, accel_beta])
        accel_fixed_frame = rotation_matrix @ accel_mobile_frame
    
    else:

        accel_r = -g0
        accel_beta = 0

        accel_mobile_frame = np.array([accel_r, accel_beta])
        accel_fixed_frame = rotation_matrix @ accel_mobile_frame

    return accel_fixed_frame

def compute_accel_fixed(beta, current_stage, pitch_angle):
    """
    Computes the acceleration of the rocket in the fixed reference frame (x,y)
    Important note : In this mode, pitch angle is the computed between the rocket thrust vector and the horizontal axis in the fixed reference frame

    Inputs :
        - beta : angle between mobile frame of reference origin and fixed reference frame origin (in radians)
        - current_stage : current stage of the rocket (0 to max_n_stage)
        - pitch_angle : pitch angle between the rocket thrust vector and the horizontal axis in the fixed reference frame (in radians)
    Outputs :
        - accel_fixed_frame : vector acceleration of the rocket in the fixed reference frame (x,y) (in m/s^2)
    """

    # Numpy rotation matrix from mobile frame to fixed frame with vector_in_fixed_frame = rotation_matrix @ vector_in_mobile_frame
    rotation_matrix_mobile_to_fixed = np.array([[np.cos(beta), -np.sin(beta)],
                                                [np.sin(beta), np.cos(beta)]])
    gravity_accel = rotation_matrix_mobile_to_fixed @ np.array([g0,0])         # In the fixed reference frame (x,y)


    if prop_mass[current_stage]>=0:
        # In the fixed reference frame

        # Acceleration due to thrust and gravity only (no air drag)
        accel_x = thrust[current_stage]/np.sum(total_mass) * np.cos(pitch_angle) - gravity_accel[0]
        accel_y = thrust[current_stage]/np.sum(total_mass) * np.sin(pitch_angle) - gravity_accel[1]

    
    else:
        # In the fixed reference frame

        accel_x = - gravity_accel[0]
        accel_y = - gravity_accel[1]

    accel_fixed_frame = np.array([accel_x, accel_y])

    return accel_fixed_frame


def compute_r(x,y):
    """Computes the distance from the center of the earth to the rocket
    
    Inputs :
        - x : x coordinate of the rocket in the fixed reference frame (in m)
        - y : y coordinate of the rocket in the fixed reference frame (in m)
    Outputs :
        - r : distance from the center of the earth to the rocket (in m)    
    """
    return np.sqrt(x**2+y**2)

def compute_alt(r):
    """Computes the altitude of the rocket
    
    Inputs :
        - r : distance from the center of the earth to the rocket (in m)
    Outputs :
        - alt : altitude of the rocket (in m)
    """
    return r - earth_radius

def compute_beta(x,y):
    """Computes the angle between the rocket and the center of the earth
    
    Inputs :
        - x : x coordinate of the rocket in the fixed reference frame (in m)
        - y : y coordinate of the rocket in the fixed reference frame (in m)
    Outputs :
        - beta : angle between the rocket and the horizontal axis in the fixed reference frame (in radians)
    """
    return np.arctan2(y,x)

def compute_speed(v, a, dt):
    """Computes the speed of the rocket relative to the center of the earth
    
    Inputs :
        - v : speed of the rocket in the fixed reference frame (in m/s)
        - a : acceleration of the rocket fixed reference frame (in m/s^2)
        - dt : time step (in s)
    Outputs :
        - v : speed of the rocket fixed reference frame (in m/s)
    """
    return v + a*dt

def compute_pos(pos, v, dt):
    """Computes the position of the rocket from the center of the earth
    
    Inputs :
        - pos : position of the rocket in the fixed reference frame (in m)
        - v : speed of the rocket in the fixed reference frame (in m/s)
        - dt : time step (in s)
    Outputs :
        - pos : position of the rocket in the fixed reference frame (in m)
        """
    return pos + dt*v

def stage_sep():
    """Separates the current stage from the rocket by setting its mass to 0"""
    struct_mass[n_stage] = 0  # Structure mass of the current stage is set to 0 (ie the stage is considered as separated)

def compute_GT_pitch_angle(v):
    """Computes the pitch angle of the rocket for gravity turn
    
    Inputs :
        - v : speed of the rocket in the fixed reference frame (in m/s)
    Outputs :
        - pitch_angle : pitch angle between the rocket thrust vector and the horizontal axis in the fixed reference frame (in radians)"""
    return np.arctan2(v[1],v[0])


    

############ Initialisation of the simulation ##############

# Mode choice
mode = 0
print("Please chose a mode for simulation :")
print("0 : Simulation with constant pitch angle variation")
print("1 : Simulation with gravity turn only")
print("2 : Simulation with constant pitch angle")
mode = int(input("Mode : "))

# Initialisation of the variables
time_tab = np.linspace(0,end_time, int(end_time/dt))    # Discretized time array
pos_tab = []                                            # Position array of the rocket in the fixed reference frame (x,y)
v_tab = []                                              # Speed array of the rocket in the fixed reference frame (x,y) 
M_tab = []                                              # Mass array of the rocket (structure + propellant)
alt_tab = []                                            # Altitude array of the rocket
v=np.array([0,0])                                       # Current speed of the rocket in the fixed reference frame (x,y), initially set to 0
pos=np.array([0,earth_radius])                          # Current position of the rocket in the fixed reference frame (x,y), initially set to the earth surface
n_stage = 0                                             # Current stage of the rocket

############ Simulation ############

if mode == 0:

    ## Pitch angle control law ##
    ## Note : Pitch angle is computed in the mobile frame of reference (r,beta) with 0 radian = toward orthoradial direction in the trigonometric direction 
    ## And + pi/2 radian = toward radial direction in the trigonometric direction

    pitch_angle_tab1 = np.linspace(0, -np.pi/2,int(np.floor(pitch_slope_time/dt)))    # Pitch angle tab for the "slope" part of the control law
    pitch_angle_tab = np.append(pitch_angle_tab1, -np.pi/2*(np.ones(int(len(time_tab)-len(pitch_angle_tab1)))))   # Pitch angle tab for the constant pitch angle part of the control law

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
            
            alt = compute_alt(r)
            beta = compute_beta(pos[0], pos[1])
            a = compute_accel_mobile(beta, n_stage, k)
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

elif mode == 1:
    pitch_angle = np.pi/2 - init_pitch_diff      # Not 0 to initiate gravity turn
    pitch_angle_tab = []

    for k in range(len(time_tab)):
        if time_tab[k] == 0:
            pitch_angle_tab.append(pitch_angle)
            pos_tab.append(pos)
            v_tab.append(np.sqrt(v[0]**2+v[1]**2))
            M_tab.append(np.sum(total_mass))
            alt_tab.append(compute_alt(compute_r(pos[0], pos[1])))

        else:
            r = compute_r(pos[0], pos[1])
            
            if r < earth_radius:
                print("Crash")
                break
            
            alt = compute_alt(r)
            beta = compute_beta(pos[0], pos[1])
            a = compute_accel_fixed(beta, n_stage, pitch_angle)
            v = compute_speed(v, a, dt)
            pitch_angle = compute_GT_pitch_angle(v)
            pos = compute_pos(pos, v, dt)

            v_tab.append(np.sqrt(v[0]**2+v[1]**2))
            pos_tab.append(pos)
            alt_tab.append(alt)
            pitch_angle_tab.append(pitch_angle)
            prop_mass[n_stage] -= m_dot[n_stage]*dt

            if prop_mass[n_stage]<= 0:
                if n_stage<max_n_stage:
                    stage_sep()
                    n_stage += 1
                else:
                    prop_mass[n_stage] = 0
            
            total_mass = prop_mass + struct_mass
            M_tab.append(np.sum(total_mass))

elif mode == 2:
    pitch_angle_tab = constant_pitch * np.ones(len(time_tab))

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
            
            alt = compute_alt(r)
            beta = compute_beta(pos[0], pos[1])
            a = compute_accel_mobile(beta, n_stage, k)
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

else:
    raise Exception("Error : Please chose a valid mode")

############ Plotting ############

pos_tab = np.array(pos_tab)
v_tab = np.array(v_tab)


fig, axes = plt.subplots(2,2)
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

ax1.plot(pos_tab[:min(int(burn_time/dt),len(pos_tab)),0], pos_tab[:min(int(burn_time/dt),len(pos_tab)),1], label='Burn', color='r') # Plotting the burn phase
ax1.plot(pos_tab[min(int(burn_time/dt),len(pos_tab)):,0], pos_tab[min(int(burn_time/dt),len(pos_tab)):,1], label='Idle', color='g') # Plotting the idle phase
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