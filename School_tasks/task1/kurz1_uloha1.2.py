import numpy as np

#Uhlova vzdialenost
def premenaDeg(Ah, Amin, Asec, Ddeg, Dmin, Dsec):
    alpha = Ah*15 + Amin*1/4 + Asec*1/240
    alpha_rad = np.radians(alpha)
    delta = Ddeg + Dmin/60 + Dsec/3600
    delta_rad = np.radians(delta)
    return alpha_rad, delta_rad

Castor = premenaDeg(7, 34, 35.87319, 31, 53, 17.816)
AlfaCentauri = premenaDeg(14, 39, 29.71993, -60, -49, -55.999)

RACastor = Castor[0]
DECCastor = Castor[1]
RAAlfaCent = AlfaCentauri[0]
DECAlfaCent = AlfaCentauri[1]

cosD = np.sin(DECAlfaCent)*np.sin(DECCastor)+np.cos(DECAlfaCent)*np.cos(DECCastor)*np.cos(RAAlfaCent-RACastor)
D = np.arccos(cosD)
D_deg = np.degrees(D)
print("Uhlova vzdialenosti medzi hviezdami je " + str(round(D_deg, 4)) + " stupnov.")
