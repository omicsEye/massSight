# core.py

import numpy as np
import pandas as pd
import polars as pl
from sklearn.linear_model import HuberRegressor
from scipy.interpolate import interp1d
from sklearn.cluster import DBSCAN
from statsmodels.gam.api import GLMGam
from .khipu import khipu, format_khipu_output

def ppm_diff(mz1: float, mz2: float) -> float:
    return abs(mz1 - mz2) / mz1 * 1e6

def find_initial_matches(ds1_mz_khipus_pd: pd.DataFrame, ds2_mz_khipus_pd: pd.DataFrame, ppm_tolerance: float = 5) -> pd.DataFrame:
    matches = []

    for _, row1 in ds1_mz_khipus_pd.iterrows():
        mz1 = row1['mz_group']
        
        # Find matches in ds2 within 5 ppm
        matches_mask = ds2_mz_khipus_pd['mz_group'].apply(lambda x: ppm_diff(mz1, x) <= ppm_tolerance)
        matching_rows = ds2_mz_khipus_pd[matches_mask]
        
        # Add matches to list
        for _, row2 in matching_rows.iterrows():
            matches.append({
                'id1': row1['id'],
                'id2': row2['id'], 
                'mz1': row1['mz_group'],
                'mz2': row2['mz_group'],
                'ppm_diff': ppm_diff(row1['mz_group'], row2['mz_group'])
            })

    # Convert matches to DataFrame
    matches_df = pd.DataFrame(matches)

    return matches_df

def get_high_confidence_matches(matches_df: pd.DataFrame, 
                                ds1: pd.DataFrame, 
                                ds2: pd.DataFrame, 
                                rt_tolerance: float = 0.5) -> pd.DataFrame:
    high_confidence_matches = matches_df.merge(
    ds1, left_on='id1', right_on='id', how='left'
    ).merge(
        ds2, left_on='id2', right_on='id', how='left'
    ).query('abs(RT_x - RT_y) < 0.5'
    ).query('(isotope_x == isotope_y)'
    ).query('(modification_x == modification_y)'
    ).assign(delta_MZ = lambda x: x.MZ_y - x.MZ_x
    ).assign(delta_RT = lambda x: x.RT_y - x.RT_x
    )

    return high_confidence_matches

def get_huber_model(high_confidence_matches: pd.DataFrame) -> tuple[HuberRegressor, HuberRegressor]:
    # get huber model for mz and rt
    mz_huber = HuberRegressor().fit(high_confidence_matches[['MZ_x']], high_confidence_matches['delta_MZ'])
    rt_huber = HuberRegressor().fit(high_confidence_matches[['RT_x']], high_confidence_matches['delta_RT'])

    return mz_huber, rt_huber


def correct_ds2(ds2: pd.DataFrame, mz_huber: HuberRegressor, rt_huber: HuberRegressor) -> pd.DataFrame:
    mz_huber_coef = mz_huber.coef_[0]
    rt_huber_coef = rt_huber.coef_[0]
    mz_huber_intercept = mz_huber.intercept_
    rt_huber_intercept = rt_huber.intercept_
    ds2["RT_corrected"] = (ds2["RT"] - rt_huber_intercept)/(rt_huber_coef + 1)
    ds2["MZ_corrected"] = (ds2["MZ"] - mz_huber_intercept)/(mz_huber_coef + 1)
    return ds2

def get_duplicate_matches(final_matches: pd.DataFrame) -> pd.DataFrame:
    duplicate_matches = final_matches.groupby('id1').size().reset_index(name='count')
    duplicate_matches = duplicate_matches[duplicate_matches['count'] > 1]

    return duplicate_matches


def get_final_matches(final_matches: pd.DataFrame, duplicate_matches: pd.DataFrame) -> pd.DataFrame:
    # Initialize list to store best matches
    best_matches = []

    # For each compound with multiple matches
    for compound_id in duplicate_matches['id1']:
        # Get all matches for this compound
        compound_matches = final_matches[final_matches['id1'] == compound_id].copy()

        # Calculate RT and MZ differences using corrected values
        compound_matches['rt_diff'] = abs(compound_matches['RT_x'] - compound_matches['RT_corrected'])
        compound_matches['mz_diff'] = abs(compound_matches['MZ_x'] - compound_matches['MZ_corrected'])
        
        # Score each match based on RT and MZ differences
        # Lower score is better
        compound_matches['match_score'] = compound_matches['rt_diff'] + compound_matches['mz_diff']
        
        # Get the match with the lowest score
        best_match = compound_matches.loc[compound_matches['match_score'].idxmin()]
        
        # Add to best matches list
        best_matches.append({
            'id1': best_match['id1'],
            'id2': best_match['id2']
        })

    # Create DataFrame of best matches
    best_matches_df = pd.DataFrame(best_matches)

    # Combine best matches with non-duplicate matches
    final_matches = pd.concat([
        final_matches[~final_matches['id1'].isin(duplicate_matches['id1'])],
        best_matches_df
    ])

    return final_matches

def mass_combine(ds1_file: str, 
                 ds2_file: str,
                 id_col: int,
                 rt_col: int,
                 mz_col: int,
                 int_cols: tuple[int, int],
                 delimiter: str = "\t",
                 mode: str = "pos",
                 ppm_tolerance: float = 5,
                 rt_tolerance: float = 0.5):

    ds1 = khipu(ds1_file, id_col=id_col, rt_col=rt_col, mz_col=mz_col, int_cols=int_cols, delimiter=delimiter, mode=mode, mz_tolerance_ppm=ppm_tolerance, rt_tolerance=rt_tolerance)
    ds2 = khipu(ds2_file, id_col=id_col, rt_col=rt_col, mz_col=mz_col, int_cols=int_cols, delimiter=delimiter, mode=mode, mz_tolerance_ppm=ppm_tolerance, rt_tolerance=rt_tolerance)

    ds1_mz_khipus_pd = format_khipu_output(ds1)
    ds2_mz_khipus_pd = format_khipu_output(ds2)

    initial_matches = find_initial_matches(ds1_mz_khipus_pd, ds2_mz_khipus_pd, ppm_tolerance=ppm_tolerance)
    high_confidence_matches = get_high_confidence_matches(initial_matches, ds1, ds2, rt_tolerance=rt_tolerance)

    mz_huber, rt_huber = get_huber_model(high_confidence_matches)
    ds2 = correct_ds2(ds2, mz_huber, rt_huber)

    matches = initial_matches.merge(
        ds1, left_on='id1', right_on='id', how='left'
        ).merge(
            ds2, left_on='id2', right_on='id', how='left'
            ).query(f'abs(RT_x - RT_corrected) < {rt_tolerance}'
            ).query(f'abs((MZ_corrected - MZ_x)/MZ_x * 1e6) < {ppm_tolerance}'
            ).query('(isotope_x == isotope_y) or (isotope_x.isna()) or (isotope_y.isna())'
            ).query('(modification_x == modification_y) or (modification_x.isna()) or (modification_y.isna())'
            ).assign(delta_MZ = lambda x: x.MZ_corrected - x.MZ_x
            ).assign(delta_RT = lambda x: x.RT_corrected - x.RT_x
            )

    duplicate_matches = get_duplicate_matches(matches)
    final_matches = get_final_matches(matches, duplicate_matches)

    return final_matches