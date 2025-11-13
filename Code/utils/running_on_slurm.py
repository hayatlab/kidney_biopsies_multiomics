def running_on_slurm():
    import os
    try:
        from IPython import get_ipython
        interactive = get_ipython() is not None
    except Exception:
        interactive = False
    slurm = bool(
        os.environ.get("SLURM_JOB_ID")
        or os.environ.get("SLURM_JOB_NODELIST")
        or os.environ.get("SLURM_NTASKS")
        or os.environ.get("SLURM_ARRAY_TASK_ID")
    )
    return slurm or (not interactive and __name__ == "__main__")
