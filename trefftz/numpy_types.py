"""
Type aliases for NumPy-based numerical arrays.

This module defines reusable type aliases for commonly used NumPy array
types throughout the Trefftz library.

The aliases improve:

- type readability
- static analysis support
- IDE autocompletion
- documentation clarity

Aliases
-------
complex_array
    NumPy array containing complex128 data.

float_array
    NumPy array containing float64 data.

int_array
    NumPy array containing int64 data.

Notes
-----
The aliases are based on :class:`numpy.typing.NDArray` and are intended
for type annotations only. They do not enforce runtime validation of
array shapes or dtypes.
"""


from numpy import float64, complex128, int64
from numpy.typing import NDArray
from typing import TypeAlias

complex_array: TypeAlias = NDArray[complex128]
float_array: TypeAlias = NDArray[float64]
int_array: TypeAlias = NDArray[int64]
