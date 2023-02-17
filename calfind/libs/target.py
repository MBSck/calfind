import numpy as np
import astropy.units as u

from typing import Optional
from astropy.time import Time
from astropy.units import Quantity
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from datetime import datetime, timedelta

from .readout import get_delay_line_restrictions


class Target:
    """"""
    def __init__(self, name: str, observatory_name: str) -> None:
        self.name = name
        self.coordinates = SkyCoord.from_name(self.name)
        self.location = EarthLocation.of_site(observatory_name)

    def get_altitude_and_azimuth(self, time: Time):
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

    # TODO: Implement margin for extended visibility of target and calibrator
    def get_sky_motion(self, start_time: Time,
                       duration: int, margin: Optional[int] = 2) -> SkyCoord:
        """Calculates azimuth for a timeframe

        Parameters
        ----------
        start_time: Time
            The start of the observation
        duration: int
            The duration of the observation
        margin: int, optional
            A margin around the start and end of the observation

        Returns
        -------
        sky_motion: SkyCoord
            An array of the altitude and azimuth for a given duration
        """
        time = np.linspace(start_time, start_time+timedelta(minutes=duration), duration*2)
        return self.get_altitude_and_azimuth(time)

    def get_airmass(self, start_time: Time,
                    duration: int, margin: Optional[int] = 2) -> Quantity:
        """Gets the airmass over a duration"""
        return self.get_sky_motion(start_time, duration, margin).secz

    def check_observability(self, time: Time, array_configuration: str):
        azimuth_and_altitude = self.get_altitude_and_azimuth(Time(time))

    # TODO: Use SkyCoord's "positional_angle" to calculate "great_circle angle"
    # TODO: Use SkyCoord's "separation" to calculate azimuthal distance?
    def compare(self, time: Time, other):
        test = self.get_altitude_and_azimuth(time)
        azimuthal_distance = test.position_angle(other)
        angular_distance = test.separation(other)


if __name__ == "__main__":
    target = Target("HD72106B", "paranal")
    print(target.get_sky_motion(Time(datetime.utcnow()), 30))
    print(type(target.get_airmass(Time(datetime.utcnow()), 30)))

