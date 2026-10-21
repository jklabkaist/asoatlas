"""ASO amenability annotation.

    from aso_amenability import annotate_variants
    resL = annotate_variants([(chrom, pos, ref, alt, transcript_id, vepD)])

The SpliceAI model is loaded on the first access to `annotate_variants` /
`annotate_variant`, not when this package is imported.
"""

import os
import sys

__version__ = '1.0.0'

# the vendored jklib / jkspliceai keep their own `import jklib` / `import jkspliceai`
_vendor_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vendor')
if os.path.isdir(_vendor_dir) and _vendor_dir not in sys.path:
    sys.path.insert(0, _vendor_dir)

_lazyD = {'annotate_variants': 'aso_amenability.annotate',
          'annotate_variant': 'aso_amenability.annotate',
          'amenability_result': 'aso_amenability.utils',
          'coding_only_result': 'aso_amenability.utils',
          'classify_amenability': 'aso_amenability.utils',
          'amenability_headerL': 'aso_amenability.utils',
          'sai_headerL': 'aso_amenability.utils',
          'sai_cons_headerL': 'aso_amenability.utils',
          'amen_headerL': 'aso_amenability.utils'}

__all__ = list(_lazyD)


def __getattr__(name):
    if name in _lazyD:
        import importlib
        return getattr(importlib.import_module(_lazyD[name]), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
