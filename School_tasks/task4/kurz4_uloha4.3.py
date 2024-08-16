from astropy.time import Time
import astropy.units as units
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz

obstime = Time(58273.916666, format='mjd', scale='utc')
coordinates = SkyCoord(x=-3792984.942, y=-3517357.760, z=5908907.466, unit='m',
    representation_type='cartesian', frame='gcrs', obstime=obstime)

location_BA = EarthLocation('48d09m03.92sN', '17d04m12.36sE', height=300*units.m)
location_BB = EarthLocation('48d44m12.56sN', '19d08m16.30sE', height=349*units.m)
location_ZahadneMestoNaS = EarthLocation('48d44m45.72sN', '22d10m58.97sE', height=120*units.m)

coord_fk5 = coordinates.fk5
print("Object's RA: " + str(coord_fk5.ra) + ", DEC: " + str(coord_fk5.dec))

coorAltazBA = coordinates.transform_to(AltAz(obstime=obstime, location=location_BA))
print("Bratislava; altitude: " + str(coorAltazBA.alt) + ", azimuth: " + str(coorAltazBA.az))

coorAltazBB = coordinates.transform_to(AltAz(obstime=obstime, location=location_BB))
print("Banska Bystrica; altitude: " + str(coorAltazBB.alt) + ", azimuth: " + str(coorAltazBB.az))

coorAltazS = coordinates.transform_to(AltAz(obstime=obstime, location=location_ZahadneMestoNaS))
print("A mysterious town which is either Snina or Sobrance; altitude: " + str(coorAltazS.alt) +
    ", azimuth: " + str(coorAltazS.az))
