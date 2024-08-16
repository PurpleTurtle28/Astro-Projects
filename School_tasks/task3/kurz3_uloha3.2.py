import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def polynom_1(x, a, b):
    return a*x + b

def polynom_2(x, a, b, c):
    return a*x**2 + b*x + c

def polynom_3(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

def polynom_4(x, a, b, c, d, e):
    return a*x**4 + b*x**3 + c*x**2 + d*x + e

def polynom_5(x, a, b, c, d, e, f):
    return a*x**5 + b*x**4 + c*x**3 + d*x**2 + e*x + f

def polynom_6(x, a, b, c, d, e, f, g):
    return a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2 +f*x + g

def residuals(observ,calc):
    dif = list()
    for i in range(n_measurements):
        dif.append(observ[i]-calc[i])
    return dif

def rms(d):
    delta = 0
    for i in range(n_measurements):
        delta += (d[i]**2)
    return np.sqrt(delta/n_measurements)
    
n_measurements = 90
x_array = np.linspace(0, 450, n_measurements)

amplitude = np.random.normal(12, 2)
shift = np.random.normal(25, 2)

y_array = np.cos(np.radians(x_array))*amplitude + shift
y_noise = 1.5 * np.random.normal(size=x_array.size)
y_array_n = y_array + y_noise

fig = plt.figure(figsize=(12,8))

fig.text(0.5, 0.04, 'x values', ha='center')
fig.text(0.04, 0.5, 'y values', va='center', rotation='vertical')

axes = fig.subplots(3,2)

poly_array = [[1,2,3], [4,5,6]]
RMS_array = list()

for i in range(0,2):
    for j in range(0,3):
        temp = "polynom_"
        temp2 = str(poly_array[i][j])
        temp3 = temp + temp2

        popt, pcov = curve_fit(locals()[str(temp3)], np.radians(x_array), y_array_n)
        axes[j,i].plot(x_array, locals()[str(temp3)](np.radians(x_array), *popt), label = "Found Fit")
        axes[j,i].plot(x_array, y_array_n, '.', label = "Measurements")

        delta = residuals(y_array_n, (locals()[str(temp3)](np.radians(x_array), *popt)))
        RMS_array.append(rms(delta))
        
plt.figure(2)
x = [1,2,3,4,5,6]
plt.plot(x, RMS_array, '*', markersize = '8')
plt.xlabel("order of polynomial function")
plt.ylabel("RMS value")

plt.show()
