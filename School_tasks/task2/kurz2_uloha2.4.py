import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import random 

def isItALine(line, columns):
    condition = False
    
    if (len(line) == columns):
        condition = True

    return condition
    

def grafy(my_file, num_of_years, period_between_years):
    year = list()
    month = list()
    day = list()
    SMA = list()        #semi major axis
    eccentricity = list() 
    inclination = list()       
    AOP1 = list()        #argument of periapsis
    RAAN1 = list()       #right ascention of the ascending node
    perihelium = list()
    longitude = list()
    mean_anomaly = list()

    with open(my_file, 'r') as in_file:
        for line in in_file:
            lin = line.split()
            if (isItALine(lin, 11) == False):
                continue
            y, m, d, S, e, i, A, R, p, l, ma = line.split()
            year.append(float(y))
            month.append(float(m))
            day.append(float(d))
            SMA.append(float(S))
            eccentricity.append(float(e))
            inclination.append(float(i))
            AOP1.append(float(A))
            RAAN1.append(float(R))
            perihelium.append(float(p))
            longitude.append(float(l))
            mean_anomaly.append(float(ma))
    
    with open("data/elements.txt", 'a+') as in_file2:
        for i in range(0, 1827, period_between_years):
            in_file2.write(str(SMA[i]) + "\t" + str(eccentricity[i]) + "\t" + str(inclination[i]) + "\t" + str(AOP1[i]) + "\t" + str(RAAN1[i]) + '\n')

    a = list()
    ecc = list()
    incl = list() 
    AOP = list()
    RAAN = list()

    with open("data/elements.txt", 'r') as in_file2:
        for line in in_file2:
            sma, e, i, w, o, = line.split()
            a.append(float(sma))
            ecc.append(float(e))
            incl.append(np.radians(float(i)))
            AOP.append(np.radians(float(w)))
            RAAN.append(np.radians(float(o)))

    fig = plt.figure(figsize=(16,8))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("z")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("z")
    
    counter = 0
    for k in range(0,7+num_of_years):
        
        counter += 1
        eccentric_anomaly = list()
        
        for i in range(0, 360, 1):
            E = [0]
            for j in range(0, int(1e4)):
                pom = (E[j]-ecc[k]*np.sin(E[j])-np.radians(i)) / (1-ecc[k]*np.cos(E[j]))
                E.append(E[j] - pom)
                if abs(E[j]-E[j+1]) < 1e-12:     
                    break
            E = (E[-1])
            eccentric_anomaly.append(E)

        x_array = list()
        y_array = list()
        P_array = list()
        Q_array = list()
        vector_r = list()

        r_x = list()
        r_y = list()
        r_z = list()

        for l in range(0, len(eccentric_anomaly)):
            x = a[k]*(np.cos(eccentric_anomaly[l])-ecc[k])
            y = a[k]*np.sqrt(1-ecc[k]**2)*np.sin(eccentric_anomaly[l])


            P = np.array(
                [[np.cos(AOP[k])*np.cos(RAAN[k])-np.sin(AOP[k])*np.cos(incl[k])*np.sin(RAAN[k])],
                [np.cos(AOP[k])*np.sin(RAAN[k])+np.sin(AOP[k])*np.cos(incl[k])*np.cos(RAAN[k])],
                [np.sin(AOP[k])*np.sin(incl[k])]])

            Q = np.array(
                [[-np.sin(AOP[k])*np.cos(RAAN[k])-np.cos(AOP[k])*np.cos(incl[k])*np.sin(RAAN[k])],
                [-np.sin(AOP[k])*np.sin(RAAN[k])+np.cos(AOP[k])*np.cos(incl[k])*np.cos(RAAN[k])],
                [np.cos(AOP[k])*np.sin(incl[k])]])

            x_array.append(x)
            y_array.append(y)
            P_array.append(P)
            Q_array.append(Q)

            vector_r.append(x*P + y*Q)
        
        for m in range(len(vector_r)):
            r_x.append(float(vector_r[m][0]))
            r_y.append(float(vector_r[m][1]))
            r_z.append(float(vector_r[m][2]))

        r = random.random()
        
        if k < 5:
            ax1.plot(r_x, r_y, r_z, color=(r, 0, 1))
        if k == 6:
            ax1.scatter(r_x, r_y, r_z, color=(r, 0, 1))     
        elif k > 6 and k <= counter:
            ax2.scatter(r_x, r_y, r_z, color=(r, 0, 1))     

    plt.show()
    
    file = open("data/elements.txt", 'r+')
    file.truncate(264)
    file.close()

    fig.savefig("figure_uloha2.4.png")

grafy("data/01501.aei", 12, 120)
