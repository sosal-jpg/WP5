import numpy as np
from numpy import sin, cos,pi

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

    def __init__(self, Name=None, E=0, y_stress=0, alpha=1, rho=1,nu=0.5):
        self.Name = Name
        self.E = E
        self.y_stress = y_stress
        self.alpha = alpha
        self.rho = rho
        self.nu = nu

mat_att=Material('ds')
mat_fast=Material('Material')
t1=0.003
F_z=1e3
F_xy=1e3

def mass(x):
    a_,b_,t2,w_,D2=x
    e=1.5*D2
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D2
    return mat_att.rho(((a+b+t2)*w-4*pi*D2**2/4)*t2)+4*mat_fast.rho*(t1 + t2 + 0.8 * D2 * 2 * 1.5 ** 2)*(pi*D2**2/4)
def sandwich_side_failure(x):
    a_,b_,t2,w_,D2=x
    return ((F_xy/(D2*t2))**2+3*)**(1/2)
    

