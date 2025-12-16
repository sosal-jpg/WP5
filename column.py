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
def ColumnBuckling(R,t,L,E):
    A = np.pi * R**2
    I = np.pi * R**3 * t
    Sigma_cb= np.pi**2 * E * I/ (A * L**2)
    return Sigma_cb


def ShellBuckling(R,t1,p,L,E,nu):
    Q= p*R**2/(E*t1**2)
    lamda = np.sqrt((12* L**4 * (1-nu**2))/(np.pi**4 * R**2 * t1**2) )
    k = lamda + 12* L**4 * (1-nu**2) / ( np.pi**4 * R**2 * t1**2 * lamda)
    Sigma_sb= (1.983 - 0.983*np.exp(-23.14*Q))*k*np.pi**2*E*t1**2/(12*(1-nu**2)*L**2)
    return Sigma_sb



  #output è total or cricical stress  
  #maximum stress che deriva da carichi effettivi quindi loads e local pressure (stress_max)
  #poi funzione per buckling compute the critical stress for this; poi fai sigma_max/ sigma_critical di buckling, ovviamente se questa 
  #spits out value maggiore di 1 allroa fails by buckling :
  # quindi due funzioni per buckling una per sheet e una per 

