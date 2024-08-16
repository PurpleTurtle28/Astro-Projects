import astropy.units as units
from astropy import constants as const

AU = const.au
dNeptune = 30.6 * AU
dPluto = 34.9 * AU
dVoyager1 = 148.6 * AU
dVoyager2 = 123.6 * AU
dAlphaCentauri = 4.4 * 63241.07709 * AU

light_sp = const.c
t_dne = dNeptune / light_sp /60 /units.s *units.min
t_dpl = dPluto / light_sp /60 /units.s *units.min
t_dv1 = dVoyager1 / light_sp /60 /units.s *units.min
t_dv2 = dVoyager2 / light_sp /60 /units.s *units.min
t_dac = dAlphaCentauri / light_sp / 86400 / 365.25 / units.s * units.year

print("Time that light needs to travel from: \nNeptune: " + str(t_dne.round(4)) + "\nPluto: " +
    str(t_dpl.round(4)) + "\nVoyager 1: " + str(t_dv1.round(4)) + "\nVoyager 2: " +
    str(t_dv2.round(4)) + "\nA. Centauri: " + str(t_dac.round(4)))