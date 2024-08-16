import numpy as np
import astropy.units as units
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, Angle
from astropy.coordinates import get_body, EarthLocation, AltAz, SkyCoord
import matplotlib.pyplot as plt
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs72
from sgp4.propagation import sgp4

def iss(tle1, tle2, mjd):

    satellite = twoline2rv(tle1, tle2, wgs72)
    jd_sat_ref = satellite.jdsatepoch
    mjd_sat_ref = jd_sat_ref - 2400000.5

    az = list()
    alt = list()

    for i in range(0, 4*60):
        delta_time = ((mjd_1+i/(24*3600)) - mjd_sat_ref) * 1440
        vec_teme_pos, vec_teme_vel = sgp4(satellite, delta_time)

        obstime = Time(mjd_1, format='mjd', scale='utc')

        coordinates = SkyCoord(x=vec_teme_pos[0], y=vec_teme_pos[1], z=vec_teme_pos[2], 
            unit='km', representation_type='cartesian', frame='gcrs', obstime=obstime)

        satellite_altaz = coordinates.transform_to(AltAz(obstime=obstime, location=my_location))
        az.append(satellite_altaz.az/units.deg)
        alt.append(satellite_altaz.alt/units.deg)

    ax1.plot(np.radians(az), alt, '--', color='salmon', label='ISS')


def nightsky(data):
    ras = list()
    decs = list()

    with open(data, 'r') as my_file:
        for line in my_file:
            ra, dec = line.split()
            ras.append(Angle(ra))
            decs.append(Angle(dec))
    
    alts = list()
    azs = list()
    for i in range(0, len(decs)):
        coord = SkyCoord(ra=ras[i], dec=decs[i], representation_type='spherical', frame='gcrs', obstime=time_1)
        vec_pos_altaz = coord.transform_to(AltAz(obstime=time_1, location=my_location))
        alts.append(vec_pos_altaz.alt)
        azs.append(vec_pos_altaz.az)

    for j in range(0, len(alts)):
        ax1.plot(np.radians(azs[j]), alts[j], color='lemonchiffon', marker='.')


if __name__ == "__main__":
    my_location = EarthLocation(lat=48.37253*units.deg, lon=17.27363*units.deg, height=536.1*units.m)

    time_1 = Time("2020-03-23 18:02:00", format='iso', scale='utc')
    time_2 = Time("2020-03-23 18:06:00", format='iso', scale='utc')
    mjd_1 = 58931.75138889
    string = str(time_1)
    string2 = str(time_2)

    tle_line1 = "1 25544U 98067A   20076.84952521  .00000812  00000-0  22770-4 0  9994"
    tle_line2 = "2 25544  51.6436  85.3376 0006075  32.6186  28.6959 15.49249183217705"

    plt.figure(figsize=(8, 8))

    ax1 = plt.subplot(polar=True)
    ax1.set_theta_zero_location('N')
    ax1.set_theta_direction(-1)
    ax1.set_rlim(90, -10)
    ax1.set_yticks([0, 90])
    ax1.set_facecolor("navy")
    ax1.set_title("Tracking ISS, AGO, " + string[0:-4] + "-" + string2[11:19])

    nightsky('stars.dat')
    iss(tle_line1, tle_line2, mjd_1)

    plt.legend()
    plt.show()
