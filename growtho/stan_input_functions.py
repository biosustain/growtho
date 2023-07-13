"""Functions for generating input to Stan from prepared data."""


from typing import Dict

import numpy as np

from growtho.data_preparation import PreparedData


def get_stan_input(ppd: PreparedData) -> Dict:
    """General function for creating a Stan input."""
    times = ppd.times
    species = ppd.coords["species"] 
    substrates = ppd.coords["substrate"] 
    species_ix = dict(zip(species, np.arange(len(species)) + 1))
    return {
        "N_species": len(ppd.coords["species"]),
        "N_reactor": len(ppd.coords["reactor"]),
        "N_timepoint": len(ppd.coords["time"]),
        "N_conc_measurement": len(ppd.conc_measurements),
        "N_biomass_measurement": len(ppd.biomass_measurements),
        "is_substrate": [1 if s in substrates else 0 for s in species],
        "reactor_yconc": ppd.conc_measurements["reactor"],
        "species_yconc": ppd.conc_measurements["species"].map(species_ix),
        "timepoint_yconc": ppd.conc_measurements["time_ix"],
        "reactor_ybiomass": ppd.biomass_measurements["reactor"],
        "timepoint_ybiomass": ppd.biomass_measurements["time_ix"],
        "timepoint_time": times["time"],
        "ybiomass": ppd.biomass_measurements["y"],
        "yconc": ppd.conc_measurements["y"]
    }

