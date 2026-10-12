"""Domain selector: PB_DOMAIN=v1 (default) or v2 picks scripts/domains/<name>.py, and every path
(inputs, results, runs) is per domain.  Everything that names a cell lives in the spec.
"""
import importlib
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAME = os.environ.get("PB_DOMAIN", "v1")
_spec = importlib.import_module(f"domains.{NAME}")
globals().update({k: v for k, v in vars(_spec).items() if not k.startswith("__")})
INPUTS = os.path.join(ROOT, "inputs", NAME)
RESULTS = os.path.join(ROOT, "results", NAME)
RUNS = os.path.join(ROOT, "runs", NAME)
NX, NY = IJ_BOUNDS[2] - IJ_BOUNDS[0], IJ_BOUNDS[3] - IJ_BOUNDS[1]


def local(conus2_i, conus2_j):
    """conus2 cell -> local (i, j) in the box."""
    return conus2_i - IJ_BOUNDS[0], conus2_j - IJ_BOUNDS[1]
