import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def star_motion(my_file):
    time = list()
    pos_x_axis = list()
    pos_y_axis = list()
    vel_x_axis = list()
    vel_y_axis = list()
    angular_mom = list()

    with open(my_file, 'r') as in_file:
        for line in in_file:
            contr_num, init_dist, init_circ_vel, t, px, py, vx, vy, am = line.split()
            contr_num = contr_num
            init_dist = init_dist
            init_circ_vel = init_circ_vel
            time.append(float(t))
            pos_x_axis.append(float(px))
            pos_y_axis.append(float(py))
            vel_x_axis.append(float(vx))
            vel_y_axis.append(float(vy))
            angular_mom.append(float(am))
    
    fig = plt.figure(figsize=(18,7))
    ax = fig.add_subplot(1,3,1, projection='3d')
    ax.plot(pos_x_axis, pos_y_axis, time, c='g')
    ax.set_xlabel("x position [kpc]")
    ax.set_ylabel("y position [kpc]")
    ax.set_zlabel("time [yrs]")
    
    total_velocity = list()
    for i in range(0, len(time)):
        total_velocity.append(np.sqrt(vel_x_axis[i]**2 + vel_y_axis[i]**2))
    
    ax = fig.add_subplot(1,3,2, projection='3d')
    ax.plot(pos_x_axis, pos_y_axis, total_velocity, c='g')
    ax.set_xlabel("x position [kpc]")
    ax.set_ylabel("y position [kpc]")
    ax.set_zlabel("total velocity [km/s]")
    
    distance = list()
    for i in range(0, len(time)):
        distance.append(np.sqrt(pos_x_axis[i]**2 + pos_y_axis[i]**2))
    
    ax = fig.add_subplot(1,3,3, projection='3d')
    ax.scatter(distance, time, total_velocity, c=total_velocity, cmap='viridis')
    ax.set_xlabel("distance [kpc]")
    ax.set_ylabel("time [yrs]")
    ax.set_zlabel("total velocity [km/s]")

    plt.show()
    fig.savefig("figure_uloha2.3.png")

star_motion("data/data5.0188a0.9557.dat")