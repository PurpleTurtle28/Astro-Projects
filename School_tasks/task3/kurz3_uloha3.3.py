import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate
import matplotlib.pyplot as plt

def residuals(observ, calc, phaseArg):
    dif = list()
    for i in range(0, len(phaseArg)):
        dif.append(observ[i]-calc[i])
    return dif

def rms(d, phaseArg):
    delta = 0
    for i in range(0, len(phaseArg)):
        delta += d[i]**2
    return np.sqrt(delta/len(phaseArg))

def fourier4(x, w, a0, a1, b1, a2, b2, a3, b3, a4, b4):    
    y = a0 + a1 * np.cos(x * w) + b1 * np.sin(x * w) + \
        a2 * np.cos(2 * x * w) + b2 * np.sin(2 * x * w) + \
        a3 * np.cos(3 * x * w) + b3 * np.sin(3 * x * w) + \
        a4 * np.cos(4 * x * w) + b4 * np.sin(4 * x * w) 
    return y

def fourier6(x, w, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6):    
    y = a0 + a1 * np.cos(x * w) + b1 * np.sin(x * w) + \
        a2 * np.cos(2 * x * w) + b2 * np.sin(2 * x * w) + \
        a3 * np.cos(3 * x * w) + b3 * np.sin(3 * x * w) + \
        a4 * np.cos(4 * x * w) + b4 * np.sin(4 * x * w) + \
        a5 * np.cos(5 * x * w) + b5 * np.sin(5 * x * w) + \
        a6 * np.cos(6 * x * w) + b6 * np.sin(6 * x * w)
    return y

def fourier8(x, w, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8):    
    y = a0 + a1 * np.cos(x * w) + b1 * np.sin(x * w) + \
        a2 * np.cos(2 * x * w) + b2 * np.sin(2 * x * w) + \
        a3 * np.cos(3 * x * w) + b3 * np.sin(3 * x * w) + \
        a4 * np.cos(4 * x * w) + b4 * np.sin(4 * x * w) + \
        a5 * np.cos(5 * x * w) + b5 * np.sin(5 * x * w) + \
        a6 * np.cos(6 * x * w) + b6 * np.sin(6 * x * w) + \
        a7 * np.cos(7 * x * w) + b7 * np.sin(7 * x * w) + \
        a8 * np.cos(8 * x * w) + b8 * np.sin(8 * x * w)
    return y

def fourier10(x, w, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10):    
    y = a0 + a1 * np.cos(x * w) + b1 * np.sin(x * w) + \
        a2 * np.cos(2 * x * w) + b2 * np.sin(2 * x * w) + \
        a3 * np.cos(3 * x * w) + b3 * np.sin(3 * x * w) + \
        a4 * np.cos(4 * x * w) + b4 * np.sin(4 * x * w) + \
        a5 * np.cos(5 * x * w) + b5 * np.sin(5 * x * w) + \
        a6 * np.cos(6 * x * w) + b6 * np.sin(6 * x * w) + \
        a7 * np.cos(7 * x * w) + b7 * np.sin(7 * x * w) + \
        a8 * np.cos(8 * x * w) + b8 * np.sin(8 * x * w) + \
        a9 * np.cos(9 * x * w) + b9 * np.sin(9 * x * w) + \
        a10 * np.cos(10 * x * w) + b10 * np.sin(10 * x * w)
    return y

def fourier12(x, w, a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a6, b6, a7, b7, a8, b8, a9, b9, a10, b10, a11, b11, a12, b12):    
    y = a0 + a1 * np.cos(x * w) + b1 * np.sin(x * w) + \
        a2 * np.cos(2 * x * w) + b2 * np.sin(2 * x * w) + \
        a3 * np.cos(3 * x * w) + b3 * np.sin(3 * x * w) + \
        a4 * np.cos(4 * x * w) + b4 * np.sin(4 * x * w) + \
        a5 * np.cos(5 * x * w) + b5 * np.sin(5 * x * w) + \
        a6 * np.cos(6 * x * w) + b6 * np.sin(6 * x * w) + \
        a7 * np.cos(7 * x * w) + b7 * np.sin(7 * x * w) + \
        a8 * np.cos(8 * x * w) + b8 * np.sin(8 * x * w) + \
        a9 * np.cos(9 * x * w) + b9 * np.sin(9 * x * w) + \
        a10 * np.cos(10 * x * w) + b10 * np.sin(10 * x * w) + \
        a11 * np.cos(11 * x * w) + b11 * np.sin(11 * x * w) + \
        a12 * np.cos(12 * x * w) + b12 * np.sin(12 * x * w)
    return y


def fun(file_name, n_fig, fit_function):

    phase = list()
    magnitude = list()
    with open(file_name, "r") as my_file:
        for line in my_file:
            if line.startswith("#"):
                pass
            else:
                ph, mag, nothing = line.split()
                phase.append(float(ph))
                magnitude.append(float(mag))

    magnitude_np = np.asarray(magnitude)
    phase_np = np.asarray(phase)

    temp = str(fit_function)

    init_c = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]

    if temp[17] == "6":
        pom = [1, 0, 1, 0]
        init_c += pom
    if temp[17] == "8":
        pom = [1, 0, 1, 0, 1, 0, 1, 0]
        init_c += pom
    if temp[17:19] == "10":
        pom = [1, 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2]
        init_c += pom
    if temp[17:19] == "12":
        pom = [1, 2, 1, 2, 1, 0, 1, 2, 1, 2, 1, 0, 1, 2, 1, 2]
        init_c += pom
    
    popt, pcov = curve_fit(fit_function, phase_np, magnitude_np, p0=init_c, maxfev=20000)

    fig = plt.figure(n_fig, figsize = (8,6))

    fig.text(0.5, 0.04, 'phase', ha='center')

    plt.subplot(211)

    plt.plot(phase_np, fit_function(phase_np, *popt), color = 'Blue', label = "Found Fit")
    plt.plot(phase_np, magnitude_np, '.', color = 'SteelBlue', label = "Measurements")
    plt.ylabel("magnitude")
    plt.legend()

    if temp[17] == "6" or "8" or "4":
        plt.title("Fit function: fourier of " + str(temp[17]) + "th order")
    if temp[17:19] == "10" or "12":
        plt.title("Fit function: fourier of " + str(temp[17:19]) + "th order")
    
    delta = residuals(magnitude_np, fit_function(phase_np, *popt), phase)

    plt.subplot(212)
    
    plt.plot(phase_np, delta, color = 'ForestGreen', label = "Residuals; RMS: " + str(round(rms(delta, phase), 4)) + " mag")
    plt.ylabel("O-C")
    
    plt.legend()
    

if __name__ == "__main__":
    fun("FoldedCurves/08022B_R_2_3_4_20180321_PHASEdata.txt", 1, fourier4)
    fun("FoldedCurves/08022B_R_2_3_4_20180321_PHASEdata.txt", 2, fourier6)
    fun("FoldedCurves/08022B_R_2_3_4_20180321_PHASEdata.txt", 3, fourier8)
    fun("FoldedCurves/08022B_R_2_3_4_20180321_PHASEdata.txt", 4, fourier10)
    fun("FoldedCurves/08022B_R_2_3_4_20180321_PHASEdata.txt", 5, fourier12)

#dalsie pouzite data - treba si nastavit cislo figure (2. argument) 
# a fitovaciu funkciu (3. argument)

    #fun("FoldedCurves/17048B_R_2_3_4_5_20180321_PHASEdata.txt", 2, fourier4)
    #fun("FoldedCurves/74039C_20171015_R_2_3_4_PHASEdata.txt", 3, fourier4)

    plt.show()
   