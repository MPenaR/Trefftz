'''Types hints for numeric arrays.'''

from numpy import floating, complexfloating, integer
from numpy.typing import NDArray
from typing import TypeAlias

complex_array: TypeAlias = NDArray[complexfloating]
float_array: TypeAlias = NDArray[floating]
int_array: TypeAlias = NDArray[integer]
