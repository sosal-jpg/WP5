class Material:
    """
    Data class to store material properties required for structural analysis.
    Attributes:
        Name (str): Material identifier.
        E (float): Young's Modulus [Pa].
        y_stress (float): Yield Strength [Pa].
        alpha (float): Thermal Expansion Coefficient [1/K].
        rho (float): Density [kg/m^3].
    """

    def __init__(self, Name=None, E=0, y_stress=0, alpha=1, rho=1):
        self.Name = Name
        self.E = E
        self.y_stress = y_stress
        self.alpha = alpha
        self.rho = rho
        
import numpy as np
from numpy import sin, cos