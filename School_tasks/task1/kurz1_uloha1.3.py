import numpy as np

#Polohovy vektor
a = 11.747335 #[AU]
e = 0.92813
i = 64.850138 #[°]
V_omega = 32.702746 #[°]
m_omega = 131.821109 #[°]
E = 88.15022 #[°]
M_S = 1.989e30 #[kg]
G = 6.67408e-11 #[m^3 kg^-1 s^-2]
AU_to_km = 149597871 #[km]

x = a*(np.cos(np.radians(E))-e)
y = a*np.sqrt(1-e**2)*np.sin(np.radians(E))

n = np.sqrt(G*M_S/(a*AU_to_km*1000)**3)
E_der = n/(1-e*np.cos(np.radians(E)))

x_der = -a*AU_to_km*E_der*np.sin(np.radians(E))
y_der = a*AU_to_km*E_der*np.sqrt(1-e**2)*np.cos(np.radians(E))
print("x = " + str(x) + " AU\ny = " + str(y) + " AU")
print("dx/dt = " + str(x_der) + " km/s\ndy/dt = " + str(y_der) + " km/s")
