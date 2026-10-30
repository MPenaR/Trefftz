class TrialFunction:
    def __mul__(self, other: "TestFunction") -> "SesquilinearForm":



class TestFunction:
    
    def __rmul__(self, other: "TrialFunction") -> "SesquilinearForm":



class SesquilinearForm:
    pass


