class Material():
    pass

class Dielectric(Material):
    def __init__(self, relative_permittivity: float):
        self.eps_r = relative_permittivity

class Metallic(Material):
    pass