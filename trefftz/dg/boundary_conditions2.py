from typing import Protocol, Optional, Any


class BoundaryCondition(Protocol):
    data: Optional[Any]


class NeumannBC:
    def __init__(self, data: Any | None = None):
        self.data = data


# class SoundHardBC:
#     def __init__(self, data: object | None = None):
#         self.data = data


class SoundHardBC:
    def __init__(self):
        self.data = None


class RadiatingBC:
    def __init__(self, data: object | None = None):
        self.data = data


class NtDBC:
    def __init__(self, truncating_radius: float, data: object | None = None):
        self.truncating_radius = truncating_radius
        self.data = data


class CircularNtD:
    ...


class WaveguideNtD:
    ...
