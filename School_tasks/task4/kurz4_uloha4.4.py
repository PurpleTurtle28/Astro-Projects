import numpy as np
from astropy.time import Time
import astropy.units as units
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_moon, get_sun
import astropy.coordinates as coord

obstime = Time('2021-04-01T0:00:00.000', format='isot', scale='utc')
location_BA = EarthLocation('48d09m03.92sN', '17d04m12.36sE', height=300*units.m)


def phases(time, location, num_of_days):
    potential_time = list()
    with open("moonphases.txt", "w") as MyData:
        MyData.write("Date\t\t\t\t\t\tMoon phase")

        i = 0
        while i < num_of_days:
            moon = coord.get_moon(time=time+i, location=location, ephemeris='builtin')    
            sun = get_sun(time + i)
            elong = sun.separation(moon)
            temp1 = sun.distance*np.sin(elong)
            temp2 = moon.distance - sun.distance*np.cos(elong)
            temp3 = np.arctan2(temp1, temp2)
            phase = (1 + np.cos(temp3))/2

            if phase > 0.995:                           # 1 = full moon, 0 = new moon
                MyData.write("\n" + str((time+i).iso) + "\t\t" + str(phase)[0:5] + "\t\tclose to full moon")
                potential_time.append(time+i)
            else:
                MyData.write("\n" + str((time+i).iso) + "\t\t" + str(phase)[0:5])
            i += 1

    eclipses(potential_time)


def eclipses(potential_time):
    for i in potential_time:
        
        j = 0
        while j < 1440:
            t = i-0.5+j/1440
            gs = get_sun(t)
            gm = get_moon(t)
            dif = abs(gs.ra - gm.ra)
            temp = str(dif)

            if temp[0:6] == '180d00':
                print(t)
                gsdec = gs.dec
                gmdec = gm.dec
                pom = abs(abs(gsdec) - abs(gmdec))/2
                
                if pom <= coord.angles.Angle('0.52d'):
                    print("Probable lunar eclipse")
                if pom >= coord.angles.Angle('0.52d'):
                    print("No lunar eclipse")
            j += 1

phases(obstime, location_BA, 90)
