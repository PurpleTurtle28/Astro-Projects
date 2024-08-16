import numpy as np

#Rotacia vektorov
i = 64.850138 #[°]
V_omega = 32.702746 #[°]
m_omega = 131.821109 #[°]

V_o = np.radians(V_omega)
m_o = np.radians(m_omega)
i_rad = np.radians(i)

P = np.array(
    [[np.cos(m_o)*np.cos(V_o)-np.sin(m_o)*np.cos(i_rad)*np.sin(V_o)],
    [np.cos(m_o)*np.sin(V_o)+np.sin(m_o)*np.cos(i_rad)*np.cos(V_o)],
    [np.sin(m_o)*np.sin(i_rad)]])

Q = np.array(
    [[-np.sin(m_o)*np.cos(V_o)-np.cos(m_o)*np.cos(i_rad)*np.sin(V_o)],
    [-np.sin(m_o)*np.sin(V_o)+np.cos(m_o)*np.cos(i_rad)*np.cos(V_o)],
    [np.cos(m_o)*np.sin(i_rad)]])

#udaje z ulohy 1.3:
x = -10.523860123661745
y = 4.370739864296843
x_der = -8.954974564322603
y_der = 0.10766008241477583

vektor_r = x*P + y*Q
vektor_der_r = x_der*P + y_der*Q
print("Zlozky r su: " + str(vektor_r[0]) + ", " + str(vektor_r[1]) + ", " + str(vektor_r[2]))
print("Zlozky dr/dt su: " + str(vektor_der_r[0]) + ", " + str(vektor_der_r[1]) + ", " + str(vektor_der_r[2]))
