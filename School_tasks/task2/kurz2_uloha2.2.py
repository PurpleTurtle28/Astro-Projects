import matplotlib.pyplot as plt
import numpy as np

def grafy(my_file, period):
    time = list()
    rel_int = list()
    int_err = list()
    magn = list()
    magn_err = list()

    with open(my_file, 'r') as in_file:
        for line in in_file:
            if line.startswith("#"):
                pass
            elif len(line.strip()) == 0:
                pass
            else:
                t, useless, ri, ie, m, me = line.split()
                useless = useless
                time.append(float(t))
                rel_int.append(float(ri))
                int_err.append(float(ie))
                magn.append(float(m))
                magn_err.append(float(me))
    
    plt.figure(figsize=(12,7))
    plt.subplot(211)
    plt.errorbar(time, magn, yerr = magn_err, fmt = 'o')
    plt.xlim([0,500])
    plt.ylim(top=17)
    plt.xlabel("time [s]")
    plt.ylabel("magnitude")

    phase = list()
    for i in time:
        temp = i/period
        phase.append(temp)
    
    modif_data = list()
    for i in phase:
        if i <= 1:
            modif_data.append(i)
        else:
            while i >1:
                i -=1
            modif_data.append(i)
    
    plt.subplot(212)
    plt.errorbar(modif_data, magn, yerr = magn_err, fmt = 'o')
    plt.xlim([0,1])
    plt.xlabel("phase")
    plt.ylabel("magnitude")
    
    plt.show()

grafy("data/14085A_2_3_5_6_7_8_9_R_20180530_DATA.txt", 110.76)