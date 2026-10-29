"""merged_timecourse — merge N developmental time-course datasets onto a shared
pseudo-time axis and derive per-gene kinetics (t50, fold change, direction,
program categories). A Python recreation of the R merged_transcript pipeline,
generalized to any number of datasets and to incremental addition.

    from merged_timecourse import Dataset, merge_datasets, add_dataset

    scott  = Dataset.from_files("scott_GNP",  "scott_normalized.tsv",  "scott_meta.tsv",  POSITIONS)
    hatten = Dataset.from_files("hatten_GNP", "hatten_normalized.tsv", "hatten_meta.tsv", POSITIONS)

    merged = merge_datasets([scott, hatten])        # job 3: N datasets
    merged = merged.add_dataset(third_dataset)      # job 2: add to an existing merge
    merged.kinetics.loc["ENSMUST00000114623"]       # t50 / fold change / program
"""
from .dataset import Dataset
from .merge import MergedTimecourse, merge_datasets, add_dataset
from . import normalize, core, plot

__all__ = [
    "Dataset", "MergedTimecourse", "merge_datasets", "add_dataset",
    "normalize", "core", "plot",
]
