from astropy.time import Time
import astropy.units as units
import astropy.coordinates as coord
from astropy.coordinates import EarthLocation

def Alt_Az(time, location, day_separation):

    with open("sunrisets.txt", "w") as MyData:
        MyData.write("Date\t\t\t\t\t\tSun's altitude")
        
        d = 0
        while d < 366:
            month = time + d
            i = 0
            while i < 850:
                time_m = month + i/(24*45)
                altaz = coord.AltAz(location=location, obstime=time_m)
                sun = coord.get_sun(time_m)
                position = sun.transform_to(altaz)
                i += 1
                gg = str(position.alt)
                if gg[0:3] == "0d0" or gg[0:4] == "-0d0":
                    MyData.write("\n" + str(time_m.iso) + "\t\t" + str(position.alt))
                
            d += day_separation

obstime = Time('2020-01-01T0:00:00.000', format='isot', scale='utc')
my_location = EarthLocation(lat=48.37253*units.deg, lon=17.27363*units.deg, height=536.1*units.m)

Alt_Az(obstime, my_location, 31)
