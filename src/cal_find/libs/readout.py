import numpy as np

from pathlib import Path


DELAY_LINES = {"small": "A0-B2-C1-D0", "medium": "K0-G2-J3-D0",
               "large": "A0-G1-J3-J2", "astrometric": "A0-G1-K0-J2",
               "UTs": "U1-U2-U3-U4"}


# NOTE: These only work for the 'paranal' location and VLTI
def get_delay_line_restrictions(array_configuration: str):
    """Gets the restrictions for the VLTI/Paranal in azimuth and altitude corresponding to the
    array_configuration

    Parameters
    ----------
    array_configuration: str
        The array configuration. Can be "small", "medium", "large", "astrometric" and
        "UTs2"

    Returns
    -------
    altitude: np.ndarray
    azimuth: np.ndarray
    """
    delay_restriction_file = Path(__file__).parents[3] / "data/delay_line_restrictions" \
            / "combined" / f"{DELAY_LINES[array_configuration]}.npy"
    return np.load(delay_restriction_file, allow_pickle=True)[::-1]


if __name__ == "__main__":
    alt, az = get_delay_line_restrictions("small")
    print(alt, az, sep="\n\n")

