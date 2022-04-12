"""Top-level imports.

See the following SO post for a disambiguation of what's going on here.

    https://stackoverflow.com/a/44842

The same pattern is used in PyMC3 for top-level `pm.name_of_things`.
"""
from .constructor import *
from .sequences import *
from .validation import *


# The functions exposed from submodules of `utils` available in `seqlike.utils`.
__all__ = (
    add_seqnums_to_letter_annotations,
    ungap,
    slice_seqrec,
)
