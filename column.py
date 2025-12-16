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
from numpy import sin, cos,pi
from scipy.optimize import minimize, NonlinearConstraint
g=9.8

bounds = [
    (0.00, 100.00),  # R
    (0.00, 0.3),  # t1
    (1.01, 5.0000),  # L
]
options = {"maxiter": 100, }
mats_col=[]
mats_fast=[]

sandwich=[[]]
axial_a=6*g
transterse_a=2*g
p=5e5
#devo iterare per il caso in cui sommo pressione e senza
for mat_col in mats_col:
    rho_col=mat_col.rho
    def choose_p(func, p):
        modified=lambda x: func(x,p)
        return modified
    def MaxStress_compressive(x,p):
        R,t1,L=x
        A=2*pi*R*t1
        I=pi*R**3*t1
        M_bend = mass(x)*transterse_a*L/2
        F_axial=mass(x)*axial_a-p*pi*R**2
        stress_hoop=p*R/t1
        for mass,h in sandwich:
            M_bend+=mass*transterse_a*h
            F_axial+=mass*axial_a
        theta=np.arange(0,91,1)
        def stress(theta):
            Q=R**2*t1*sin(theta)
            stress_z=F_axial/A+M_bend*R/I
            stress_shear=F_axial*Q/(I*t1)
            return (stress_z**2+stress_hoop**2-stress_z*stress_hoop+3*stress_shear**2)**(1/2)
        return max(stress(i) for i in theta)

    def ColumnBuckling(x,p):
        R,t1,L=x 
        A = 2 * pi * R * t1
        I = pi * R**3 * t1
        Sigma_cb= pi**2 * mat_col.E * I/ (A * L**2)
        return MaxStress_compressive(x)/Sigma_cb # devo fare divisione con maximal stress nel grafico e devo controllare che questo non superi 1, penso che questo maximal stress derivi da loading 
    def ShellBuckling(x,p):
        R,t1,L=x
        Q= p*R**2/(mat_col.E*t1**2)
        lamda = np.sqrt((12* L**4 * (1-mat_col.nu**2))/(pi**4 * R**2 * t1**2) )
        k = lamda + 12* L**4 * (1-mat_col.nu**2) / ( pi**4 * R**2 * t1**2 * lamda)
        Sigma_sb= (1.983 - 0.983*np.exp(-23.14*Q))*k*pi**2*mat_col.E*t1**2/(12*(1-mat_col.nu**2)*L**2)
        return MaxStress_compressive(x)/Sigma_sb
    def mass(x):
        R,t1,L=x
        return rho_col*(pi*R*t1*2*L)
    
    constraints=[]
    constraints.append(NonlinearConstraint(choose_p(ColumnBuckling,0), lb=1,ub=np.inf))
    constraints.append(NonlinearConstraint(choose_p(ShellBuckling,0), lb=1,ub=np.inf))
    constraints.append(NonlinearConstraint(choose_p(MaxStress_compressive,0), lb=0,ub=mat_col.y_stress))
    constraints.append(NonlinearConstraint(choose_p(ColumnBuckling,p), lb=1,ub=np.inf))
    constraints.append(NonlinearConstraint(choose_p(ShellBuckling,p), lb=1,ub=np.inf))
    constraints.append(NonlinearConstraint(choose_p(MaxStress_compressive,p), lb=0,ub=mat_col.y_stress))
    x0=[1.2,0.003,3]
    res=minimize(mass, x0=x0, constraints=constraints, bounds=bounds, method='SLSQP', options=options)



  #output è total or cricical stress  
  #maximum stress che deriva da carichi effettivi quindi loads e local pressure (stress_max)
  #poi funzione per buckling compute the critical stress for this; poi fai sigma_max/ sigma_critical di buckling, ovviamente se questa 
  #spits out value maggiore di 1 allroa fails by buckling :
  # quindi due funzioni per buckling una per sheet e una per 


