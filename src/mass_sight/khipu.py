import pandas as pd
import numpy as np
import contextlib
import io
from khipu.extended import *

def khipu(
    file_path: str,
    id_col: int,
    rt_col: int,
    mz_col: int,
    int_cols: tuple[int, int],
    delimiter: str = "\t",
    mode: str = "pos",
    mz_tolerance_ppm: float = 5,
    rt_tolerance: float = 2,
    initial_isotopes: list = [],
    initial_adducts: list = [],
    extended_adducts: list = [],
) -> pd.DataFrame:
    """
    Load data from a file using khipu, run khipu, and return a pandas DataFrame.
    """

    if initial_isotopes == []:
        # Isotope patterns for non-tracing HILIC data
        initial_isotopes = isotope_search_patterns[:3]

    if initial_adducts == []:
        # Conservative adduct patterns for initial network building
        initial_adducts = [
        # (21.9820, 'Na/H'),   # [M+Na]+
        (37.955882, 'K/H'),    # [M+K]+
        ]

    if extended_adducts == []:
        extended_adducts = [
        (17.02655, 'NH3'),
        (-18.0106, '-H2O'),
        (-17.02655, '-NH3'),
        (-18.0106, '-H2O'),
        (46.00548, 'CO2H2'),
        ]

    full_df = pd.read_csv(file_path, delimiter=delimiter)
    
    with contextlib.redirect_stdout(io.StringIO()):
        peaklist = read_features_from_text(
            open(file_path).read(),
            id_col=id_col,
            mz_col=mz_col,
            rtime_col=rt_col,
            intensity_cols=int_cols,
            delimiter=delimiter,
        )

    with contextlib.redirect_stdout(io.StringIO()):
        khipu_list, _ = peaklist_to_khipu_list(
            peaklist,
            isotope_search_patterns=initial_isotopes,
            adduct_search_patterns=initial_adducts,
            extended_adducts=extended_adducts,
            mz_tolerance_ppm=mz_tolerance_ppm,
            rt_tolerance=rt_tolerance,
            mode=mode,
            charges=[1],
        )

    kpdict, fdict = {}, {}
    for KP in khipu_list:
        kpdict[KP.id] = KP
        for f in KP.nodes_to_use:
            fdict[f] = KP

    df = {}

    for id, model in kpdict.items():
        for metabolite in model.format_to_epds()["MS1_pseudo_Spectra"]:
            df[metabolite["id"]] = {
                "Compound_ID": metabolite["id"],
                "isotope": metabolite["isotope"],
                "modification": metabolite["modification"],
                "khipu_id": id,
            }

    df = pd.DataFrame(df.values())
    full_df = full_df.merge(df, on="Compound_ID", how="left")

    return full_df

def format_khipu_output(khipu_df: pd.DataFrame) -> pd.DataFrame:
    # Split khipu_id column and handle empty values
    khipu_df['id'] = khipu_df['khipu_id'].str.split('_').str[0]
    khipu_df['mz_group'] = [float(i) if pd.notna(i) else khipu_df.iloc[idx]['MZ'] - 1.007825 for idx, i in enumerate(khipu_df['khipu_id'].str.split('_').str[1])]

    # for each NaN id, assign it a unique incremented id
    # get max id from both datasets (ids are of form kp123)
    all_valid_ids = pd.concat([khipu_df['id'].dropna(), khipu_df['id'].dropna()])
    # Filter out any non-string values and extract numeric parts
    valid_kp_ids = all_valid_ids[all_valid_ids.astype(str).str.startswith('kp')]
    if len(valid_kp_ids) > 0:
        max_id = int(valid_kp_ids.str[2:].astype(int).max())
    else:
        max_id = 0  # Default if no valid kp IDs found

    # assign unique NaN ids with incremented values
    khipu_df_na_mask = khipu_df['id'].isna()

    # For ds1, assign unique incremented IDs starting from max_id + 1
    khipu_df_na_count = khipu_df_na_mask.sum()
    if khipu_df_na_count > 0:
        new_ids_khipu_df = [f"kp{max_id + i + 1}" for i in range(khipu_df_na_count)]
        khipu_df.loc[khipu_df_na_mask, 'id'] = new_ids_khipu_df

    return khipu_df

