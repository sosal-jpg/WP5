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

L_tanks=1.8
m_tanks=210
L=L_tanks
R=0.6/2
bounds = [
    (1, 30),  # t1
]
options = {"maxiter": 1000, }
mats_col = [
    Material(Name='4130 Steel', E=205e9, y_stress=725e6, alpha=12.6e-6, rho=7850),
    Material(Name='8630 Steel', E=200e9, y_stress=550e6, alpha=11.7e-6, rho=7850),
    Material(Name='2014-T6 Aluminium', E=73.1e9, y_stress=414e6, alpha=24.4e-6, rho=2800),
    Material(Name='2024-T3 Aluminium', E=73.1e9, y_stress=345e6, alpha=24.7e-6, rho=2780),
    Material(Name='2024-T4 Aluminium', E=73.1e9, y_stress=345e6, alpha=24.7e-6, rho=2780),
    Material(Name='AZ91C-T6 Magnesium', E=44.8e9, y_stress=145e6, alpha=26e-6, rho=1810),
    Material(Name='7075-T6 Aluminium', E=73.1e9, y_stress=503e6, alpha=23.6e-6, rho=2810),
    Material(Name='2219-T851 Aluminium', E=73.1e9, y_stress=352e6, alpha=22.3e-6, rho=2840),
    Material(Name='Al-Li 2195', E=74e9, y_stress=448e6, alpha=23e-6, rho=2710),
    Material(Name='Ti-6Al-4V', E=119e9, y_stress=1200e6, alpha=9.1e-6, rho=4512),
    Material(Name='8552 Carbon-Epoxy (Hexcel)', E=190e9, y_stress=3310e6, alpha=6e-6, rho=1300),
]


mats_fast=[]

sandwich=[[6.27+2,0],[37+2,0.6],[7.77+2,1.2],[2,1.8]]
axial_a=6*g
transterse_a=2*g
p_given=5e5
#devo iterare per il caso in cui sommo pressione e senza

mass_min=np.inf
t1_min=None
mat_min=None
p_min=None
for mat_col in mats_col:
    for p in [0,p_given]:
        t1=0.0002
        def choose_p(func, p):
            modified=lambda x: func(x,p)
            return modified
        def MaxStress(t1,p):
            A=2*pi*R*t1
            I=pi*R**3*t1
            M_bend = mass(t1)*transterse_a*L/2
            F_axial=mass(t1)*axial_a-p*pi*R**2+m_tanks*axial_a
            stress_hoop=p*R/t1
            for m,h in sandwich:
                M_bend+=m*transterse_a*h
                F_axial+=m*axial_a
            M_bend+=m_tanks*transterse_a*L_tanks/2
            stress_z=-F_axial/A
            return max(((stress_z+M_bend*R/I)**2+stress_hoop**2-(stress_z+M_bend*R/I)*stress_hoop)**(1/2),((stress_z-M_bend*R/I)**2+stress_hoop**2-(stress_z-M_bend*R/I)*stress_hoop)**(1/2))/mat_col.y_stress
            
        
        def MaxStress_compr(t1,p):
            A=2*pi*R*t1
            I=pi*R**3*t1
            M_bend = mass(t1)*transterse_a*L/2
            F_axial=mass(t1)*axial_a-p*pi*R**2+m_tanks*axial_a
            for m,h in sandwich:
                M_bend+=m*transterse_a*h
                F_axial+=m*axial_a
            M_bend+=m_tanks*transterse_a*L_tanks/2
            stress_z=max(F_axial/A+M_bend*R/I,0)
            return stress_z

        def ColumnBuckling(t1,p):
            A = 2 * pi * R * t1
            I = pi * R**3 * t1
            Sigma_cb= pi**2 * mat_col.E * I/ (A * L**2)
            return MaxStress_compr(t1,p)/Sigma_cb# devo fare divisione con maximal stress nel grafico e devo controllare che questo non superi 1, penso che questo maximal stress derivi da loading 
        def ShellBuckling(t1,p):
            Q= p*R**2/(mat_col.E*t1**2)
            lamda = np.sqrt((12* L**4 * (1-mat_col.nu**2))/(pi**4 * R**2 * t1**2) )
            k = lamda + 12* L**4 * (1-mat_col.nu**2) / ( pi**4 * R**2 * t1**2 * lamda)
            Sigma_sb= (1.983 - 0.983*np.exp(-23.14*Q))*k*pi**2*mat_col.E*t1**2/(12*(1-mat_col.nu**2)*L**2)
            return MaxStress_compr(t1,p)/Sigma_sb
        def mass(t1):
            return mat_col.rho*(2*pi*R*t1*L)
    
        while True:
            if  choose_p(ColumnBuckling,p)(t1)<=1 and choose_p(ShellBuckling,p)(t1)<=1 and choose_p(MaxStress,p)(t1)<=1:
                break
            t1+=0.0001
        
        if mass(t1)<mass_min:
            mass_min=mass(t1)
            mat_min=mat_col
            p_min=p
            t1_min=t1

t1=t1_min
mat_col=mat_min
p=p_min
print(f"{mat_col.Name}: {mass_min}")
print(f't= {t1}')
print(f'pressure: {p}')
print(choose_p(ColumnBuckling,p)(t1))
print(choose_p(ShellBuckling,p)(t1))
print(choose_p(MaxStress,p)(t1))

  #output è total or cricical stress  
  #maximum stress che deriva da carichi effettivi quindi loads e local pressure (stress_max)
  #poi funzione per buckling compute the critical stress for this; poi fai sigma_max/ sigma_critical di buckling, ovviamente se questa 
  #spits out value maggiore di 1 allroa fails by buckling :
  # quindi due funzioni per buckling una per sheet e una per 


