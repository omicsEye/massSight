import pandas as pd
import numpy as np
import contextlib
import io
from khipu.extended import *

def khipu(file_path: str, 
          id_col: int,
          rt_col: int,
          mz_col: int,
          int_cols: tuple[int, int],
          delimiter: str = "\t",
          mode: str = "pos",
          mz_tolerance_ppm: float = 5,
          rt_tolerance: float = 2) -> pd.DataFrame:
    """
    Load data from a file using khipu, run khipu, and return a pandas DataFrame.
    """

    full_df = pd.read_csv(file_path, delimiter=delimiter)
    
    peaklist = read_features_from_text(
        open(file_path).read(),
        id_col=id_col, 
        mz_col=mz_col, 
        rtime_col=rt_col, 
        intensity_cols=int_cols, 
        delimiter=delimiter
    )
    
    with contextlib.redirect_stdout(io.StringIO()):
        khipu_list, all_assigned_peaks = peaklist_to_khipu_list(peaklist, 
                    isotope_search_patterns=isotope_search_patterns, 
                    adduct_search_patterns=adduct_search_patterns,
                    mz_tolerance_ppm=mz_tolerance_ppm,
                    rt_tolerance=rt_tolerance, 
                    mode=mode)
    
    kpdict, fdict = {}, {}
    for KP in khipu_list:
        kpdict[KP.id] = KP
        for f in KP.nodes_to_use:
            fdict[f] = KP
    
    df = {}

    for id, model in kpdict.items():
        for metabolite in model.format_to_epds()['MS1_pseudo_Spectra']:
            df[metabolite['id']] = {
                'Compound_ID': metabolite['id'],
                'isotope': metabolite['isotope'],
                'modification': metabolite['modification'],
                'khipu_id': id,
            }

    df = pd.DataFrame(df.values())
    full_df = full_df.merge(df, on="Compound_ID", how="left")

    return full_df
