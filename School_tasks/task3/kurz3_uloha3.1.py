import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func_cos_1(x,a):
    return a*np.cos(x)

def func_cos_2(x,d):
    return np.cos(x) + d

def func_cos_3(x, a, d):
    return a * np.cos(x) + d
    
def residuals(observ,calc):
    dif = list()
    for i in range(n_measurements):
        dif.append(observ[i]-calc[i])
    return dif

def rms(d):
    delta = 0
    for i in range(n_measurements):
        delta += d[i]**2
    return np.sqrt(delta/n_measurements)

n_measurements = 110
x_array = np.linspace(0, 450, n_measurements)

amplitude = np.random.normal(12, 2)
shift = np.random.normal(25, 2)

y_array = amplitude * np.cos(np.radians(x_array)) + shift
y_noise = 1.5 * np.random.normal(size=x_array.size)
y_array_n = y_array + y_noise

def function(cosine_func, figure_num):
    plt.figure(figure_num)
    plt.subplot(211)

    popt, pcov = curve_fit(cosine_func, np.radians(x_array), y_array_n)
    plt.plot(x_array, cosine_func(np.radians(x_array), *popt), label = "Found Fit")
    plt.plot(x_array, y_array_n, '.', label = "Measurements")
    plt.legend()

    delta = residuals(y_array_n, cosine_func(np.radians(x_array), *popt))
    print("Total RMS for " + str(figure_num) + ". function: " + str(rms(delta)))
    
    plt.subplot(212)
    plt.plot(x_array, delta, 'g--')
    plt.xlabel("x values")
    plt.ylabel("O-C residuals")

function(func_cos_1, 1)
function(func_cos_2, 2)
function(func_cos_3, 3)

plt.show()