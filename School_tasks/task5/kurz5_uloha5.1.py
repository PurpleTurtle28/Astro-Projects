import numpy as np
import jplephem, random
import astropy.units as units
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, Angle
from astropy.coordinates import get_body, EarthLocation, AltAz, SkyCoord
import matplotlib.pyplot as plt


def planets(time, location):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.gca(projection='polar')
    ax.set_theta_zero_location('W')
    ax.set_theta_direction(-1)

    for i in range(0, len(bodies)):
        name = bodies[i]
        temp = get_body(bodies[i], time, location)
        coord = temp.transform_to(AltAz(obstime=time, location=location))

        if i == 0:
            ax.plot(np.radians(coord.az-90*units.deg), coord.alt, color='yellow', marker='o', markersize = 10, label=name)
        else:
            ax.plot(np.radians(coord.az-90*units.deg), coord.alt, marker='o', color=colors[i-1], label=name)
    
    ax.set_title("Positions of planets and the Sun from AGO, " + str(time.iso))
    plt.legend(loc='lower left')
    plt.show()


def options(axes, hemisphere):
    axes.set_theta_zero_location('N')
    axes.set_yticklabels([])
    axes.set_facecolor("navy")
    axes.grid(False)

    if hemisphere == "N":
        axes.set_theta_direction(-1)
        axes.set_title("Northern hemisphere")
        axes.set_rlim(90, 0, 10)
        axes.set_yticks(np.linspace(90, 0, 10))
    elif hemisphere == "S":
        axes.set_theta_direction(1)
        axes.set_title("Southern hemisphere")
        axes.set_rlim(-90, 0, -10)
        axes.set_yticks(np.linspace(-90, 0, 10))


def nightsky(data, time, location):
    ras = list()
    decs = list()

    with open(data, 'r') as my_file:
        for line in my_file:
            ra, dec = line.split()
            ras.append(Angle(ra))
            decs.append(Angle(dec))
    
    plt.figure(figsize=(13, 7))
    
    ax1 = plt.subplot(121, polar=True)
    ax2 = plt.subplot(122, polar=True)
    options(ax1, 'N')
    options(ax2, 'S')

    for i in range(0, len(ras)):
        if decs[i] < 0:
            ax2.plot(np.radians(ras[i]), decs[i], color='lemonchiffon', marker='.')
        else:
            ax1.plot(np.radians(ras[i]), decs[i], color='lemonchiffon', marker='.')

    for i in range(0, len(bodies)):
        temp = get_body(bodies[i], time, location)
        if temp.dec < 0:
            ax2.plot(np.radians(temp.ra), temp.dec, color=colors[i-1], marker='o')
        else:
            ax1.plot(np.radians(temp.ra), temp.dec, color=colors[i-1], marker='o')


    plt.show()


if __name__ == "__main__":
    obstime = Time("2020-03-23 18:00", format='iso', scale='utc')
    my_location = EarthLocation(lat=48.37253*units.deg, lon=17.27363*units.deg, height=536.1*units.m)

    bodies = solar_system_ephemeris.bodies[1:5] + solar_system_ephemeris.bodies[6:]
    colors = ['lightgray', 'dimgray', 'peru', 'red', 'orange', 'khaki', 'turquoise', 'blue']

    planets(obstime, my_location)
    nightsky('stars.dat', obstime, my_location)
