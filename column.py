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
        
import numpy as np
from numpy import sin, cos
from scipy.optimize import minimize, NonlinearConstraint

bounds = [
    (0.00, 100.00),  # R
    (0.00, 0.3),  # t1
    (1.01, 5.0000),  # L
]
options = {"maxiter": 100, }
mats_col=[]
mats_fast=[]
for mat_col in mats_col:
    rho_col=mat_col.rho
    def ColumnBuckling(x):
        R,t1,L=x
        A = 2 * np.pi * R * t1
        I = np.pi * R**3 * t1
        Sigma_cb= np.pi**2 * mat_col.E * I/ (A * L**2)
        return Sigma_cb
    def ShellBuckling(x):
        R,t1,L=x
        Q= p*R**2/(mat_col.E*t1**2)
        lamda = np.sqrt((12* L**4 * (1-mat_col.nu**2))/(np.pi**4 * R**2 * t1**2) )
        k = lamda + 12* L**4 * (1-mat_col.nu**2) / ( np.pi**4 * R**2 * t1**2 * lamda)
        Sigma_sb= (1.983 - 0.983*np.exp(-23.14*Q))*k*np.pi**2*mat_col.E*t1**2/(12*(1-mat_col.nu**2)*L**2)
        return Sigma_sb
    def mass(x):
        R,t1,L=x
        return rho_col*(np.pi*R*t1*2*L)
    constraints=[]
    constraints.append(NonlinearConstraint(ColumnBuckling, lb=0,ub=mat_col.y_stress))
    x0=[1.2,0.003,3]
    minimize(mass, x0=x0, constraints=constraints, bounds=bounds, method='SLSQP', options=options)



  #output è total or cricical stress  
  #maximum stress che deriva da carichi effettivi quindi loads e local pressure (stress_max)
  #poi funzione per buckling compute the critical stress for this; poi fai sigma_max/ sigma_critical di buckling, ovviamente se questa 
  #spits out value maggiore di 1 allroa fails by buckling :
  # quindi due funzioni per buckling una per sheet e una per 


