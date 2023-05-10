import numpy as np
import astropy.units as u

from scipy.special import j1
from datetime import datetime, timedelta
from typing import List, Optional
from astropy.time import Time
from astropy import coordinates
from astropy.coordinates import SkyCoord, EarthLocation

import datetime
from collections import namedtuple
import ephem
import pytz
from astropy.time import Time
from astropy.coordinates import EarthLocation


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


def calculate_twilight(date: Union[str, datetime],
                       site: Optional[str] = "paranal",
                       latitude: Optional[float] = None,
                       longitude: Optional[float] = None,
                       twilight_kind: Optional[str] = "astronomical"
                       ) -> Tuple[namedtuple, namedtuple]:
    """Calculates the sunset and sunrise for the input type of twilight.

    Parameters
    ----------
    date : str or datetime.date or datetime.datetime
        Can either be a string of the form "yyyy-mm-dd" or
        a `datetime.date`/`datetime.datetime`.
    site : str, optional
        The site at which the sunset and sunrise are to be
        determined.
    latitude : float, optional
        The latitude at which the sunset and sunrise are to be
        determined [deg].
    longitude : float, optional
        The longitude at which the sunset and sunrise are to be
        determined [deg].
    twilight_kind : str, optional
        The type of twilight. Can be "civil" for -6 degrees, "nautical"
        for -12 degrees and astronomical for -18 degrees.

    Returns
    -------
    sunset : namedtuple of str
        The sunset at the twilight specified. This is a named tuple
        containing string entries for "utc", "lst" and "cet".
    sunrise : namedtuple of str
        The sunrise at the twilight specified. This is a named tuple
        containing string entries for "utc", "lst" and "cet".
    """
    TimesTuple = namedtuple("TimesTuple", ["utc", "lst", "cet"])
    if isinstance(date, (datetime.date, datetime.datetime)):
        date = date.strftime('%Y-%m-%d')

    if latitude is not None:
        location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg)
    else:
        location = EarthLocation.of_site(site)
    observer = ephem.Observer()
    observer.date = date
    observer.lat = str(location.lat.value)
    observer.lon = str(location.lon.value)

    if twilight_kind == "civil":
        observer.horizon = "-6"
    elif twilight_kind == "nautical":
        observer.horizon = "-12"
    elif twilight_kind == "astronomical":
        observer.horizon = "-18"

    sunset = observer.next_setting(ephem.Sun(), use_center=True)
    sunrise = observer.next_rising(ephem.Sun(), use_center=True)

    sunset, sunrise = map(lambda x: Time(x.datetime()), [sunset, sunrise])
    sunset_utc, sunrise_utc = map(lambda x: x.utc, [sunset, sunrise])
    sunrise_utc += datetime.timedelta(days=1)

    utc_timezone = pytz.timezone("UTC")
    sunset_aware_utc = utc_timezone.localize(sunset_utc.to_datetime())
    sunrise_aware_utc = utc_timezone.localize(sunrise_utc.to_datetime())

    sunset_lst = sunset.sidereal_time("mean", longitude=location.lon)
    sunrise_lst = sunrise.sidereal_time("mean", longitude=location.lon)

    cet_timezone = pytz.timezone("CET")
    sunset_cet = sunset_aware_utc.astimezone(cet_timezone)
    sunrise_cet = sunrise_aware_utc.astimezone(cet_timezone)

    sunset = TimesTuple(str(sunset_utc),
                        str(sunset_lst),
                        f"{sunset_cet.date()} {sunset_cet.time()}")
    sunrise = TimesTuple(str(sunrise_utc),
                         str(sunrise_lst),
                         f"{sunrise_cet.date()} {sunrise_cet.time()}")
    return sunset, sunrise


# TODO: Reimplement this as classes to take into account all three types utc, lst and
# so at once
def calculate_night_lengths(date: Union[str, datetime],
                            observation_slot: str,
                            site: Optional[str] = "paranal",
                            latitude: Optional[float] = None,
                            longitude: Optional[float] = None,
                            twilight_kind: Optional[str] = "astronomical"
                            ) -> namedtuple:
    """Calculates the length of the total night as well as the 

    Parameters
    ----------
    date : str or datetime.date or datetime.datetime
        Can either be a string of the form "yyyy-mm-dd" or
        a `datetime.date`/`datetime.datetime`.
    observation_slot : str
        This is the time that is alloted for the observations
        on the night in question. Can be any combination of the
        following "0.6h1", "1.2h2" or "1.0n". The floats can be
        any number and the string specifies the start "h1" for
        the beginning of the night, "h2" for the second half of 
        the night and "n" for the full night.
    site : str, optional
        The site at which the sunset and sunrise are to be
        determined.
    latitude : float, optional
        The latitude at which the sunset and sunrise are to be
        determined [deg].
    longitude : float, optional
        The longitude at which the sunset and sunrise are to be
        determined [deg].
    twilight_kind : str, optional
        The type of twilight. Can be "civil" for -6 degrees, "nautical"
        for -12 degrees and astronomical for -18 degrees.

    Returns
    -------
    night_duration : namedtuple
    """
    NightDuration = namedtuple("NightDuration", ["duration", "start", "end",
                                                 "observation"])
    Observation = namedtuple("Observation", ["type", "start", "end",
                                             "duration", "minutes"])
    sunset, sunrise = calculate_twilight(date, site, latitude,
                                         longitude, twilight_kind)
    sunset_utc = datetime.datetime.strptime(sunset.utc, "%Y-%m-%d %H:%M:%S.%f")
    sunrise_utc = datetime.datetime.strptime(sunrise.utc, "%Y-%m-%d %H:%M:%S.%f")
    total_duration = sunrise_utc - sunset_utc

    pattern = r"(\d+\.\d+|\d+\.\d*|\d+)([a-z]+)(\d*)"
    match = re.match(pattern, observation_slot)
    multiplicative_factor = float(match.group(1))
    part_identifier = match.group(2)
    half_night_identifier = int(match.group(3)) if match.group(3) else None

    if part_identifier == "n":
        observation_duration = total_duration
        observation_start = sunrise
    elif part_identifier == "h":
        observation_duration = (total_duration / 2) * multiplicative_factor
        observation_minutes = observation_duration.seconds // 60
        if half_night_identifier == 1:
            observation_start = sunset
            observation_end = sunset_utc + observation_duration
        elif half_night_identifier == 2:
            observation_start = sunrise_utc - observation_duration
            observation_end = sunrise_utc
    observation = Observation(observation_slot, str(observation_start),
                              str(observation_end), str(observation_duration),
                              observation_minutes)
    return NightDuration(str(total_duration), sunset.lst, sunrise.lst, observation)


if __name__ == "__main__":
    date = datetime.date(2023, 5, 4)
    sunset, sunrise = calculate_twilight(date)
    breakpoint()


if __name__ == "__main__":
    ...

