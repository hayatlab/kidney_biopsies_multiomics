import numpy as np
import pybiomart
import matplotlib.pyplot as plt


def labels_on_cluster(adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None):

    from adjustText import adjust_text
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}
    medians = {}
    adata = adata[~adata.obs.index.duplicated(keep='first')]
    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)
    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs)
                 for k, (x, y) in medians.items()]
    return adjust_text(texts, **adjust_kwargs)
