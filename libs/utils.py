import numpy as np
import astropy.units as u

from typing import List, Optional
from scipy.special import j1
from datetime import date, datetime, timedelta
from astropy.time import Time
from astropy import coordinates
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


# TODO: Maybe change this class to accept more than one target?
class Target:
    """"""
    def __init__(self, name: str, observatory_name: str) -> None:
        self.name = name
        self.coordinates = SkyCoord.from_name(self.name)
        self.location = EarthLocation.of_site(observatory_name)

    def get_altitude_and_azimuth(self, time: str):
        """

        Parameters
        ----------
        time: str

        Returns
        -------
        azimuth_and_altitude: SkyCoord[u.deg, u.deg]
        """
        return self.coordinates.transform_to(AltAz(obstime=Time(time),
                                                   location=self.location))

    # TODO: Think of passing the whole target for the hours and such? Not only the airmass?
    # TODO: Check that the midnight does not go back by an entire day, but is always the
    # coming one
    # TODO: Find a way not to use time here for midnight
    # TODO: Add timespan here
    def get_time_interval(self, period: Optional[List[int]] = [-4, 4]):
        """Calculates azimuth for a timeframe

        Parameters
        ----------
        period: List[int], optional
            The time period where the target azimuth should be calculated

        Returns
        -------
        night: SkyCoord
            The coordinates of the target for a night
        """
        sampling = np.diff(period)[0]//2*100 if np.diff(period)[0]//2*100 >= 100 else 100
        return Time(get_midnight())+np.linspace(*period, sampling)*u.hour

    def check_observability(self, time, delay_restrictions):
        azimuth_and_altitude = self.get_altitude_and_azimuth(time)
        ...

    def get_airmass(self, period: Optional[List[int]] = [-4, 4]):
        return self.get_night(period).sezc

    def compare(self, other: Target):
        ...


@u.quantity_input
def _convert_azimuth(azimuth: u.deg) -> u.deg:
    """Converts the azimuth from west-of-south to east-of-north or vice-versa

    Parameters
    ----------
    azimuth: u.deg

    Returns
    -------
    azimuth: u.deg
    """
    returns (azimuth + 180*u.deg) % 360*u.deg


def get_target_right_ascension_and_declination(target_name: str,
                                               observatory_name: str,
                                               time: str) -> SkyCoord:
    """
    Convert Horizon (Alt-Az) coordinates to Hour Angle and Declination.

    Parameters
    ----------
    target_name: str
        The name of the target as given in Simbad
    observatory_name: str
        The name of the observatory. For full list see 'EarthLocation.get_site_names()'
        Examples: "paranal", "lasilla", "Mt Graham", "Observatorio de Calar Alto", etc.
    time: str

    Returns
    -------
    azimuth_and_altitude: SkyCoord[u.deg, u.deg]
    """
    ...

def convert_utc_to_julian_datetime(time: str) -> np.float64:
    """Convert UTC-datetime to Julian-datetime

    Parameters
    ----------
    time: str

    Returns
    -------
    julian_datetime: np.float64
    """
    return Time(time, format="isot", scale="utc").jd


# TODO: Think of how to implement the time in the functions
def convert_utc_to_mean_sidereal_time(observatory_name: str,
                                      time: str) -> coordinates.angles.Longitude:
    """Converts UTC-datetime to the mean sidereal time at the given location

    Parameters
    ----------
    observatory_name: str
        The name of the observatory. For full list see 'EarthLocation.get_site_names()'
        Examples: "paranal", "lasilla", "Mt Graham", "Observatorio de Calar Alto", etc.
    time: str

    Returns
    -------
    mean_sidereal_time: coordinates.angles.Longitude
    """
    location = EarthLocation.of_site(observatory_name)
    return Time(time, format="isot", scale="utc", location=location).sidereal_time("mean")


@u.quantity_input
def calibrator_visibility(wavelength: u.um,
                          baseline: u.m, diameter: u.mas) -> u.dimensionless_unscaled:
    """Calculates the visibility of a calibrator on a given baseline,
    assuming the calibrator is a uniform disk of the given diameter,
    for each point in the wavelength grid.

    Parameters
    ----------
    wavelength: u.um
    baseline: u.m
    diameter: u.mas

    Returns
    -------
    calibrator_visibility: u.dimensionless_unscaled
    """
    product = diameter.to(u.rad).value*np.pi*baseline/wavelength.to(u.m)
    return 2*j1(product)/product


def wrap_mean(azimuth: List | np.ndarray) -> float:
    """Gives proper mean azimuth, also if source passes through North
    (and az jumps from 0 to 360)

    Parameters
    ----------
    azimuth: List | np.ndarray

    Returns
    -------
    mean_azimuth: np.ndarray
    """
    dx = (azimuth-np.roll(azimuth, 1))[1:-1]
    if np.where(dx > 300)[0].size != 0:
        azimuth[azimuth > 180] -= 360
        mean_azimuth = np.mean(azimuth)
        mean_azimuth[mean_azimuth < 0.] += 360
        return mean_azimuth
    else:
        return np.mean(azimuth)


def get_midnight():
    """Gets the midnight of UTC. Either the coming midnight if after 12pm or the midnight
    of else the midnight ahead"""
    if datetime.utcnow().hour < 12:
        return datetime.combine(datetime.utcnow(), datetime.min.time())
    return datetime.combine(datetime.utcnow()+timedelta(1), datetime.min.time())


if __name__ == "__main__":
    target = Target("HD72106B", "paranal")
    print(target.get_altitude_and_azimuth(get_midnight()).separation())


