import sqlite3
import numpy as np
import pandas as pd
import scipy.io as sio

from typing import List
from pathlib import Path
from datetime import datetime


def convert_delay_line_restrictions(conversion_path: Path):
    """Reads the delay line restriction from a (.hzn)-file out and saves the output as a
    (.npy)-file with a comment in form of a (.txt)-file

    Parameters
    ----------
    conversion_path: Path
        The path where the files for conversion reside
    """
    files_to_convert = conversion_path.glob("*.hzn")
    for hzn_file in files_to_convert:
        print(f"Working on file: {hzn_file}")
        with open(hzn_file, "r+") as hzn_f:
            azimuth, altitude, comments = [], [], []
            for line in hzn_f:
                if line.startswith("#"):
                    comments.append(line.split("#")[1].strip()+"\n")
                else:
                    az, alt = line.strip().split()
                    azimuth.append(az)
                    altitude.append(alt)
        azimuth, altitude = map(lambda x: np.array(x, dtype=float), (azimuth, altitude))
        with open(conversion_path / (Path(hzn_file).stem+"_info.txt"), "w+") as txt_file:
            txt_file.write("Converted (.hzn)-file to (.npy)-array storage\n")
            txt_file.write(f"Date of conversion: {datetime.utcnow()} UTC\n")
            txt_file.write("Format is Azimuthal [deg East-of-North] and Altitude [deg]\n")
            txt_file.write("Readout with 'az, alt = np.load(<file_name>)'\n")
            txt_file.write("------------------------\n")
            txt_file.write("INFO OLD\n")
            txt_file.write("------------------------\n")
            txt_file.writelines(comments)
        npy_file = conversion_path / (Path(hzn_file).stem+".npy")
        np.save(npy_file, np.array([azimuth, altitude], dtype=object))

        # Note: Test-file load
        azimuth, altitude = np.load(npy_file, allow_pickle=True)

def convert_dat_to_mongo_db(conversion_path):
    save_data = sio.readsav(conversion_path)
    df = pd.DataFrame(save_data["calcat"])
    print(df.columns)


class SQLiteDatabase:
    """A context manager for the 'sqlite3'-package that automatically closes the
    connection and returns dictionaries instead of tuples"""
    def __init__(self, database_path: Path) -> None:
        """"""
        self.database_path = database_path
        self.connection, self.row_factory, self.cursor = None, None, None

    def __enter__(self):
        """"""
        self.connection = sqlite3.connect(database_path)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback) -> None:
        """"""
        self.connection.close()

    def create_table(self, title: str, columns: List[str]) -> None:
        """Creates a table in the database

        Parameters
        ----------
        title: str
            The table's title
        columns: List[str]
            The table's columns
        """
        create_table_command = f"CREATE TABLE {title.lower()}"\
                f"({', '.join([column.lower() for column in columns])})"
        self.cursor.execute(create_table_command)


if __name__ == "__main__":
    conversion_path = Path(__file__).parents[3] \
            / "data" / "calibrator_catalogues" / "midi_cat_2019.dat"
    # convert_dat_to_mongo_db(conversion_path)
    sql = SQLiteDatabase("hello")

