from typing import List, Optional
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from .utils import get_midnight
from .readout import get_delay_line_restrictions


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

    # TODO: Use SkyCoord's "space_motion" to calculate visibilities for the night
    def get_time_interval(self, start_time: Time,
                          end_time: Time, margin: Optional = None):
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
        ...

    def check_observability(self, time, array_configuration: str):
        azimuth_and_altitude = self.get_altitude_and_azimuth(time)

    def get_airmass(self, period: Optional[List[int]] = [-4, 4]):
        ...

    # TODO: Use SkyCoord's "positional_angle" to calculate "great_circle angle"
    # TODO: Use SkyCoord's "separation" to calculate azimuthal distance?
    def compare(self, other):
        ...


if __name__ == "__main__":
    ...

