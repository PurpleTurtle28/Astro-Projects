import numpy as np
import matplotlib.pyplot as plt
from shapely import Point, Polygon

min_mass = 0.1
max_mass = 100
pc_to_meters = 3.086e16
year_to_seconds = 3.154e7
mSun_to_kilograms = 1.99e30

G_si = 6.67430e-11
G_new = G_si * pc_to_meters**3 / year_to_seconds**2 / mSun_to_kilograms 
#G_new = 4.30091e-3

def compute_a(N, M_min, M_max):
    ap = -2.35+1
    return N*ap/(M_max**ap-M_min**ap)

def IMF(m, m_zero, m_last, a):
    value = a*m**(-2.35)
    max = (a*m_zero**(-2.35))
    min = (a*m_last**(-2.35))
    return value, max, min

def basic_IMF(m, a):
    return a*m**(-2.35)

def random_points_within(poly, num_points):
    points = []
    min_x, min_y, max_x, max_y = poly.bounds
    
    while len(points) < num_points:
        random_point = Point([np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)])
        if (random_point.within(poly)) and random_point.y <= basic_IMF(random_point.x, compute_a(num_of_stars, min_mass, max_mass)):
            points.append(random_point)
          
    return points

def spatial_distribution(n):
    # Generating distances to approx. 10 pc
    rMagn = []
    positions_xyz = []
    for _ in range(n):
        q = np.random.uniform(0.007, 0.5)
        temp = 2*np.log(1/q - 1)
        rMagn.append(temp)

    for i in range(len(rMagn)):
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2*np.pi)
        
        x = (rMagn[i] * np.sin(theta) * np.cos(phi))
        y = (rMagn[i] * np.sin(theta) * np.sin(phi))
        z = (rMagn[i] * np.cos(theta))

        temp_app = [x, y, z]
        positions_xyz.append(temp_app)

    return positions_xyz, rMagn

def velocity_distribution(n):
    vMagn = []
    velocitiesV = [] 
    for _ in range(n):
        value_kms = 10 * np.random.uniform(0.98, 1.02)
        # Converting to pc/yr
        value_pcyr = value_kms * 1.022e-6
        vMagn.append(value_pcyr)

    for i in range(len(vMagn)):
        vx, vy, vz = generate_vector_components(vMagn[i])
        temp_app = [vx, vy, vz]
        velocitiesV.append(temp_app)

    return velocitiesV, vMagn

def generate_vector_components(magnitude):
    x = np.random.normal()
    y = np.random.normal()
    z = np.random.normal()
    
    # Scale the components to match the desired magnitude
    normalization_factor = magnitude / np.linalg.norm([x, y, z])
    x *= normalization_factor
    y *= normalization_factor
    z *= normalization_factor
    
    return x, y, z

def compute_motion(N, array_of_rs, array_of_vs, array_of_ms, softening, years, step):
    G = G_new 
    r = array_of_rs
    v = array_of_vs
    m = array_of_ms 
    current = 1

    print("initial virial: ")
    energy(masses[0:-1], initialV, initialR, G)

    counts = int(years/step)
    all_positions = np.zeros([N, counts, 3], dtype=float)
    all_velocities = np.zeros([N, counts, 3], dtype=float)

    counter = 0
    while current <= years:
        vels_after_each_iter = []
        poss_after_each_iter = []
        for i in range(N): 
            acc = ([0, 0, 0])
            for j in range(N):
                if i == j:
                    continue
                else:
                    dx = r[j][0] - r[i][0]
                    dy = r[j][1] - r[i][1]
                    dz = r[j][2] - r[i][2]
                    dist = (dx**2 + dy**2 + dz**2 + softening**2)**(3/2)
                    acc[0] += G * m[j] * (dx / dist)
                    acc[1] += G * m[j] * (dy / dist)
                    acc[2] += G * m[j] * (dz / dist)
        
            # ---------LEAPFROG METHOD-----------
            #temp_vel = kick(v, i, step, acc)
            #new_rs = new_positions(r, i, step, temp_vel)
            # -----------------------------------

            # -------7-STEP-LEAPFROG METHOD------
            new_rs, new_vel = leapfrog7(i, step, r, v, acc)
            # -----------------------------------          

            all_positions[i][counter] = [new_rs[i][0], new_rs[i][1], new_rs[i][2]]
            tempR = (new_rs[i][0]**2 + new_rs[i][1]**2 + new_rs[i][2]**2)**(1/2)
            poss_after_each_iter.append(tempR)

            # ---------LEAPFROG METHOD-----------
            #new_vel = kick(temp_vel, i, step, acc)
            # -----------------------------------

            all_velocities[i][counter] = [new_vel[i][0], new_vel[i][1], new_vel[i][2]] 
            tempV = (new_vel[i][0]**2 + new_vel[i][1]**2 + new_vel[i][2]**2)**(1/2)
            vels_after_each_iter.append(tempV)

        #print("all_positions: ", all_positions, "\n")
        #print("all_velocities: ", all_velocities, "\n")

        #energy(masses, vels_after_each_iter, poss_after_each_iter, G)
        
        current += step
        counter += 1
    
    return all_positions, counts

def kick(vel, i, dt, acc):
    vel[i][0] += acc[0] * dt / 2
    vel[i][1] += acc[1] * dt / 2
    vel[i][2] += acc[2] * dt / 2

    return vel

def new_positions(pos, i, dt, vel):
    pos[i][0] += vel[i][0] * dt
    pos[i][1] += vel[i][1] * dt
    pos[i][2] += vel[i][2] * dt 

    return pos

def leapfrog7(i, dt, pos, vel, acc):
    w = 2**(1/3)
    f = 2 - w

    pos[i][0] += vel[i][0] * dt/(2*f)
    pos[i][1] += vel[i][1] * dt/(2*f)
    pos[i][2] += vel[i][2] * dt/(2*f)

    vel[i][0] += acc[0] * dt/f
    vel[i][1] += acc[1] * dt/f
    vel[i][2] += acc[2] * dt/f
    
    pos[i][0] += vel[i][0] * (1-w)*dt/(2*f)
    pos[i][1] += vel[i][1] * (1-w)*dt/(2*f)
    pos[i][2] += vel[i][2] * (1-w)*dt/(2*f)

    vel[i][0] += acc[0] * -dt*w/f
    vel[i][1] += acc[1] * -dt*w/f
    vel[i][2] += acc[2] * -dt*w/f

    pos[i][0] += vel[i][0] * (1-w)*dt/(2*f)
    pos[i][1] += vel[i][1] * (1-w)*dt/(2*f)
    pos[i][2] += vel[i][2] * (1-w)*dt/(2*f)

    vel[i][0] += acc[0] * dt/f
    vel[i][1] += acc[1] * dt/f
    vel[i][2] += acc[2] * dt/f

    pos[i][0] += vel[i][0] * dt/(2*f)
    pos[i][1] += vel[i][1] * dt/(2*f)
    pos[i][2] += vel[i][2] * dt/(2*f)

    return pos, vel

def energy(m, v, r, G):
    E_k = 0
    E_p = 0
    
    for i in range(len(m)):
        E_k += 1/2 * m[i] * v[i]**2
    
    for i in range(len(m)-1):
        E_p += G * m[i] * m[i+1] / np.abs(r[i+1]-r[i])

    print("E_p =", E_p, ", E_k =", E_k)
    print("virial =", 2 * E_k/E_p, "\n")

def plottingPositions(n, n_of_plots, all_pos, lim):
    step = int(n/n_of_plots)
    for i in range(0, n, step):
        ax2 = plt.figure(figsize=(6,6)).add_subplot(projection='3d')
        ax2.set_xlabel("x [pc]")
        ax2.set_ylabel("y [pc]")
        ax2.set_zlabel("z [pc]")
        ax2.set_xlim(-lim, lim)
        ax2.set_ylim(-lim, lim)
        ax2.set_zlim(-lim, lim)

        for j in range(num_of_objects):
            if masses[j] <= 0.4:
                ax2.scatter(all_pos[j, i, 0], all_pos[j, i, 1], all_pos[j, i, 2], alpha=0.4, color="C0", s=15)
            elif masses[j] > 0.4 and masses[j] < 1 and masses[j] != masses[-1]:
                ax2.scatter(all_pos[j, i, 0], all_pos[j, i, 1], all_pos[j, i, 2], alpha=0.4, color="blue", s=20)
            elif masses[j] >= 1 and masses[j] != masses[-1]:
                ax2.scatter(all_pos[j, i, 0], all_pos[j, i, 1], all_pos[j, i, 2], alpha=0.4, color="darkblue", s=25)
            elif masses[j] == masses[-1] and masses[j] > 100:
                ax2.scatter(all_pos[j, i, 0], all_pos[j, i, 1], all_pos[j, i, 2], alpha=0.8, color="black", s=50)
        
        plt.title("simulation duration: " + str(i) + " years")
        plt.savefig(f"Nbody/no_{i}.png", dpi=200)
        plt.close()

if __name__ == "__main__":
    # Number of stars and allowed masses defined here
    num_of_stars = 1000

    # Creating the initial mass function
    masses_func = np.linspace(min_mass, max_mass, 100)
    imfs, max_imf, min_imf = IMF(masses_func, float(masses_func[0]), float(masses_func[-1]), compute_a(num_of_stars, min_mass, max_mass))

    # Generating possible masses within a polygon instead of the whole field to save some time
    poly = Polygon([(min_mass, min_imf), (min_mass, max_imf), (max_mass, min_imf)])
    points = random_points_within(poly, num_of_stars)
    xs = [point.x for point in points]
    ys = [point.y for point in points]

    # Plotting the initial mass function vs mass with generated points
    plt.figure(figsize=(6,6))
    plt.plot(masses_func, imfs)
    plt.scatter(xs, ys, alpha=0.4, label="generated mass points")
    plt.xlabel('M_sun', fontsize=12)
    plt.ylabel('IMF', fontsize=12)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(loc="upper right")

    masses = [i for i in xs]
    positions, initialR = spatial_distribution(num_of_stars)
    velocities, initialV = velocity_distribution(num_of_stars)

    # Plotting the velocity distribution vs radius of the star system
    plt.figure(figsize=(8,4))
    plt.scatter(initialR, initialV, alpha=0.7)
    plt.xlabel("initial positions [pc]", fontsize=12)
    plt.ylabel("initial velocities [pc/yr]", fontsize=12)
    plt.ylim([0.5e-5, 1.8e-5])

    # Adding a black hole to the center
    masses.append(1e3)
    #masses.append(1)
    positions.append([0,0,0])
    velocities.append([0,0,0])
    num_of_objects = num_of_stars+1

    # Plotting initial positions of all objects in 3D space 
    ax = plt.figure(figsize=(8,8)).add_subplot(projection='3d')
    ax.set_title("Initial positions")
    ax.set_xlim([-10, 10])
    ax.set_ylim([-10, 10])
    ax.set_zlim([-10, 10])

    for i in range(0, num_of_objects):
        if masses[i] <= 0.4:
            ax.scatter(positions[i][0], positions[i][1], positions[i][2], alpha=0.4, color="C0", s=15)
        elif masses[i] > 0.4 and masses[i] < 1 and masses[i] != masses[-1]:
            ax.scatter(positions[i][0], positions[i][1], positions[i][2], alpha=0.4, color="blue", s=20)
        elif masses[i] >= 1 and masses[i] != masses[-1]:
            ax.scatter(positions[i][0], positions[i][1], positions[i][2], alpha=0.4, color="darkblue", s=25)
        elif masses[i] == masses[-1] and masses[i] > 100:
            ax.scatter(positions[i][0], positions[i][1], positions[i][2], alpha=1, color="black", s=50)

    radius_of_cluster = 10 #pc
    soft = 1/10 * (radius_of_cluster/num_of_objects**(1/3))

    # Lauching the simulation
    all_positions, counts = compute_motion(num_of_objects, positions, velocities, masses, soft, 1000, 1)
    plottingPositions(counts, 100, all_positions, 10)

    plt.show()
