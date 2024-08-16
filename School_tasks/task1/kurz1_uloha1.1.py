import numpy as np

#Kvadraticka rovnica
def najdiKorene(a,b,c):
    x1 = ((-b)+np.sqrt(b**2-4*a*c))/(2*a)
    x2 = ((-b)-np.sqrt(b**2-4*a*c))/(2*a)
    return x1, x2

x1, x2 = najdiKorene(2,3,-2)
print("Korene kvadratickej rovnice su: " + str(najdiKorene(2,3,-2)))