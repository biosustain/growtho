"""Provides functions prepare_data_x.

These functions should take in a dataframe of measurements and return a
PreparedData object.

"""
import json
import os

import numpy as np
import pandas as pd
import pandera as pa
from pandera.typing import DataFrame, Series
from pandera import Field
from pydantic import BaseModel
from growtho import util

NAME_FILE = "name.txt"
COORDS_FILE = "coords.json"
CONC_FILE = "timecourse_conc.csv"
BIOMASS_FILE = "timecourse_biomass.csv"
TIMES_FILE = "times.csv"

HERE = os.path.dirname(__file__)
DATA_DIR = os.path.join(HERE, "..", "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PREPARED_DIR = os.path.join(DATA_DIR, "prepared")
RAW_DATA_FILES = {
    "raw_timecourse": os.path.join(RAW_DIR, "timecourse.csv"),
}


def prepare_data():
    """Run main function."""
    print("Reading raw data...")
    raw_data = {
        k: pd.read_csv(v, index_col=None) for k, v in RAW_DATA_FILES.items()
    }
    data_preparation_functions_to_run = [prepare_data_hooman]
    print("Preparing data...")
    for dpf in data_preparation_functions_to_run:
        print(f"Running data preparation function {dpf.__name__}...")
        prepared_data = dpf(raw_data["raw_timecourse"])
        output_dir = os.path.join(PREPARED_DIR, prepared_data.name)
        print(f"\twriting files to {output_dir}")
        if not os.path.exists(PREPARED_DIR):
            os.mkdir(PREPARED_DIR)
        write_prepared_data(prepared_data, output_dir)


class Measurements(pa.SchemaModel):
    """A PreparedData should have a dataframes like this for concs, and another
    one for biomass.

    Other columns are also allowed!
    """

    reactor: Series[int] = Field(ge=1)
    strain: Series[int] = Field(ge=1)
    species: Series[str]
    time_ix: Series[int] = Field(ge=1)
    y: Series[float]

class Times(pa.SchemaModel):
    """A PreparedData should have a time like this.

    Other columns are also allowed!
    """

    time: Series[float]
    time_ix: Series[int] = Field(ge=1)


class PreparedData(BaseModel, arbitrary_types_allowed=True):
    """What prepared data looks like in this analysis."""

    name: str
    coords: util.CoordDict
    conc_measurements: DataFrame[Measurements]
    biomass_measurements: DataFrame[Measurements]
    times: DataFrame[Times]


def load_prepared_data(directory: str) -> PreparedData:
    """Load prepared data from files in directory."""
    with open(os.path.join(directory, COORDS_FILE), "r") as f:
        coords = json.load(f)
    with open(os.path.join(directory, NAME_FILE), "r") as f:
        name = f.read()
    conc_measurements = pd.read_csv(os.path.join(directory, CONC_FILE))
    biomass_measurements = pd.read_csv(os.path.join(directory, BIOMASS_FILE))
    times = pd.read_csv(os.path.join(directory, TIMES_FILE))
    return PreparedData(
        name=name,
        coords=coords,
        conc_measurements=DataFrame[Measurements](conc_measurements),
        biomass_measurements=DataFrame[Measurements](biomass_measurements),
        times=DataFrame[Times](times)
    )


def write_prepared_data(prepped: PreparedData, directory):
    """Write prepared data files to a directory."""
    if not os.path.exists(directory):
        os.mkdir(directory)
        prepped.conc_measurements.to_csv(os.path.join(directory, CONC_FILE))
        prepped.biomass_measurements.to_csv(
            os.path.join(directory, BIOMASS_FILE)
        ) 
        prepped.times.to_csv(os.path.join(directory, TIMES_FILE))
    with open(os.path.join(directory, COORDS_FILE), "w") as f:
        json.dump(prepped.coords, f)
    with open(os.path.join(directory, NAME_FILE), "w") as f:
        f.write(prepped.name)

def prepare_data_hooman(timecourse_raw: pd.DataFrame) -> PreparedData:
    """Prepare data with an interaction column."""
    substrates = ["gluc"]
    products = ["lac"]
    reactor_to_strain = {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2}
    species = substrates + products
    species_cols_in = [s.capitalize() + "_adj" for s in species]
    reactor = timecourse_raw["Sample ID"].apply(
        lambda x: x.split(" ")[0].split("-")[1]
    ).astype(int).rename("reactor")
    strain = reactor.map(reactor_to_strain).rename("strain")
    biomass = 10000000 * (
        timecourse_raw["Live (cells/ml)_adj"]
        * 4/3
        * 1000
        * np.pi
        * timecourse_raw["Estimated cell diameter (um)"]/2/1000000 ** 3
    ).rename("biomass")
    time_arr = np.sort(timecourse_raw["time"].unique())
    times = pd.DataFrame({
        "time": time_arr, "time_ix":np.arange(len(time_arr))+1
    })
    timecourse_prepped = (
        timecourse_raw
        .join(reactor)
        .join(biomass)
        .join(strain)
        .join(times.set_index("time"), on="time")
        .rename(columns=dict(zip(species_cols_in, species)))
    )
    conc_measurements, biomass_measurements = (
        timecourse_prepped
        .melt(
            id_vars=["reactor", "strain", "time_ix"],
            value_vars=value_vars,
            value_name="y",
            var_name="species"
        )
        for value_vars in (species, ["biomass"])
    )
#    conc_measurements = conc_measurements.loc[lambda df: df["species"] == "lac"]
    return PreparedData(
        name="hooman",
        coords=util.CoordDict(
            {
                "reactor": reactor.astype(str).unique().tolist(),
                "strain": strain.astype(str).unique().tolist(),
                "substrate": map(lambda s: s.lower(), substrates),
                "product": map(lambda s: s.lower(), products),
                "species": map(lambda s: s.lower(), conc_measurements["species"].unique()),
                "time": times["time"].astype(str).tolist(),
            }
        ),
        conc_measurements=DataFrame[Measurements](conc_measurements),
        biomass_measurements=DataFrame[Measurements](biomass_measurements),
        times=DataFrame[Times](times)
    )
