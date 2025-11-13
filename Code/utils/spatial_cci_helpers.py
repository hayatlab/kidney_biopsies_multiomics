from pygam import LinearGAM, s
import numpy as np
import re
import os
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess


def load_continuous_chunks(prefix, folder=".", start=0, end=198361):
    """
    Load and concatenate a sequence of parquet files that form a continuous
    coverage for a given prefix.

    Files are expected to follow the naming pattern: "{prefix}_{start}_{end}.parquet".
    The function will:
      - discover files in `folder` that match the pattern,
      - sort them by their start coordinate,
      - select the files that form a continuous coverage starting at `start`
        and stopping when the coverage reaches or exceeds `end`,
      - read each parquet file into a pandas.DataFrame and concatenate them.

    Parameters
    ----------
    prefix : str
        Filename prefix to match before the coordinate range.
    folder : str, optional
        Directory to search for parquet files (default: ".").
    start : int, optional
        Expected start coordinate for coverage (default: 0).
    end : int, optional
        Target end coordinate for coverage (default: 198361).

    Returns
    -------
    pandas.DataFrame
        Concatenated DataFrame of all selected parquet files. If no files are found
        an empty DataFrame or a concatenation error from pandas may occur.
        If coverage is incomplete a warning is printed.
    """
    pattern = re.compile(rf"{prefix}_(\d+)_(\d+)\.parquet$")
    files = []

    # sort by start
    files.sort(key=lambda x: x[0])

    # ensure continuity
    continuous = []
    cur = start
    for s, e, path in files:
        if s == cur:
            continuous.append(path)
            cur = e
        if cur >= end:
            break

    if cur < end:
        print(f"Warning: coverage incomplete (stopped at {cur})")

    # load all
    dfs = [pd.read_parquet(f) for f in continuous]
    return pd.concat(dfs, ignore_index=True)


# Pivot dataframes to have gene_a as rows, gene_b as columns, score as values
def pivot_cci_score_df(df, genes_a, genes_b):
    """
    Pivot a cell–cell interaction (CCI) score DataFrame into a matrix form.

    The input DataFrame is expected to contain at least the columns:
    "gene_a", "gene_b", and "score". The function constructs a pivot table
    with rows indexed by `genes_a`, columns by `genes_b`, and values given by
    the mean of "score" where multiple entries exist. Missing combinations are
    filled with zeros and the final matrix is reindexed to include all provided
    genes in the requested order.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing "gene_a", "gene_b", and "score" columns.
    genes_a : sequence
        Sequence of gene identifiers to use as the row index (order preserved).
    genes_b : sequence
        Sequence of gene identifiers to use as the columns (order preserved).

    Returns
    -------
    pandas.DataFrame
        Pivoted DataFrame with index `genes_a` and columns `genes_b`. Missing
        values are filled with 0.
    """
    pivot = df.pivot_table(
        index="gene_a",
        columns="gene_b",
        values="score",
        aggfunc="mean",
        fill_value=0
    )
    # Ensure all genes are present in rows and columns
    pivot = pivot.reindex(index=genes_a, columns=genes_b, fill_value=0)
    return pivot


def compute_loess(x, y, frac=0.4):
    """Compute LOESS-smoothed line for sparkline.

    Parameters
    ----------
    x : array-like
        Independent variable values. Can be unsorted; lowess will sort internally.
    y : array-like
        Dependent variable values, same length as `x`.
    frac : float, optional
        The fraction of the data used when estimating each y-value (default 0.4).

    Returns
    -------
    x_s : numpy.ndarray
        Sorted and normalized x-coordinates of the smoothed curve in the range [0, 1].
    y_s : numpy.ndarray
        Corresponding smoothed y-values. If the input has fewer than 3 points,
        the raw x and y are returned unchanged.
    """
    if len(x) < 3:
        return x, y  # too few points → raw
    lo = lowess(y, x, frac=frac, return_sorted=True)
    x_s = (lo[:, 0] - lo[:, 0].min()) / (lo[:, 0].max() - lo[:, 0].min())
    y_s = lo[:, 1]
    return x_s, y_s


def compute_spline(x, y, n_splines=20, lam=0.6, return_ci=False):
    """
    Fit a P-spline with fixed smoothing (lambda).

    Parameters:
        x, y      : 1D numpy arrays of numeric data. x is normalized internally.
        n_splines : int, optional
                    Number of B-spline basis functions to use (default 20).
        lam       : float, optional
                    Smoothing penalty (lambda) applied to the GAM (default 0.6).
        return_ci : bool, optional
                    If True, return 95% confidence interval bounds along with
                    the smoothed fit (default False).

    Returns:
        If return_ci is False:
            x_norm, y_smooth
                x_norm : numpy.ndarray
                    x values normalized to the range [0, 1].
                y_smooth : numpy.ndarray
                    Predicted/smoothed y-values at x_norm.
        If return_ci is True:
            x_norm, y_smooth, ci_low, ci_high
                ci_low, ci_high : numpy.ndarray
                    Lower and upper bounds of the 95% confidence intervals.

    Notes
    -----
    The function uses pygam.LinearGAM with a single spline term. If the input
    has fewer than 3 points it returns the raw inputs (and None for CI bounds
    when requested).
    """
    if len(x) < 3:
        return (x, y) if not return_ci else (x, y, None, None)

    x = np.array(x)
    y = np.array(y)

    # Normalize x
    x_norm = (x - x.min()) / (x.max() - x.min())

    # TO REMOVE:
    # Fit GAM with fixed lambda
    # gam = LinearGAM(s(0, n_splines=n_splines)).fit(
    #     X=x_norm.reshape(-1, 1),
    #     y=y,
    #     weights=None,
    #     lam=lam
    # )
    # END TO REMOVE
    
    gam = LinearGAM(s(0, n_splines=20), lam=lam)
    gam.fit(x_norm.reshape(-1, 1), y)

    # Predict spline
    y_smooth = gam.predict(x_norm.reshape(-1, 1))

    if return_ci:
        ci = gam.confidence_intervals(x_norm.reshape(-1, 1))
        return x_norm, y_smooth, ci[:, 0], ci[:, 1]

    return x_norm, y_smooth


def compute_spline_auto(x, y, n_splines=20, return_ci=False):
    """
    Fit a P-spline with automatic smoothing (lambda) selection.

    The function normalizes x to the [0, 1] interval, fits a pygam LinearGAM
    with `n_splines` basis functions and performs a grid search to select the
    smoothing parameter automatically.

    Parameters
    ----------
    x, y : array-like
        Independent and dependent variables (1D). Converted to numpy arrays.
    n_splines : int, optional
        Number of B-spline basis functions (default: 20).
    return_ci : bool, optional
        If True, return 95% confidence interval bounds along with the smoothed fit.

    Returns
    -------
    If return_ci is False:
        x_norm, y_smooth
            x_norm : numpy.ndarray
                x values normalized to the range [0, 1].
            y_smooth : numpy.ndarray
                Predicted/smoothed y-values at x_norm.
    If return_ci is True:
        x_norm, y_smooth, ci_low, ci_high
            ci_low, ci_high : numpy.ndarray
                Lower and upper bounds of the 95% confidence intervals.
    """
    if len(x) < 3:
        return (x, y) if not return_ci else (x, y, None, None)

    x = np.array(x)
    y = np.array(y)

    # Normalize x
    x_norm = (x - x.min()) / (x.max() - x.min())

    # Fit GAM with automatic lambda search
    gam = LinearGAM(s(0, n_splines=n_splines)).gridsearch(
        X=x_norm.reshape(-1, 1),
        y=y
    )

    # Predict spline
    y_smooth = gam.predict(x_norm.reshape(-1, 1))

    if return_ci:
        ci = gam.confidence_intervals(x_norm.reshape(-1, 1))
        return x_norm, y_smooth, ci[:, 0], ci[:, 1]

    return x_norm, y_smooth
