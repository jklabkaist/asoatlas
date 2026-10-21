"""Reference files used by the amenability annotation.

Each path falls back to the copy shipped with the vendored jkspliceai / jklib
packages and can be overridden through the matching environment variable.
"""

import importlib.util
import os


def _package_dir(name):
    spec = importlib.util.find_spec(name)
    if spec is None or spec.origin is None:
        raise ImportError(f"{name} is not importable; run setup_data.sh first")
    return os.path.dirname(spec.origin)


JKSPLICEAI_DIR = _package_dir('jkspliceai')
JKLIB_DIR = _package_dir('jklib')

REFFLAT_MANE_PATH = os.environ.get(
    'ASO_AMEN_REFFLAT',
    os.path.join(JKSPLICEAI_DIR, 'annotations', 'gencode.v41.annotation.gtf.refFlat.MANE.v1.arranged.txt'))

LABRANCHOR_PATH = os.environ.get(
    'ASO_AMEN_LABRANCHOR',
    os.path.join(JKSPLICEAI_DIR, 'annotations', 'labranchor.gencode_v40.top.sorted.bed.gz'))
