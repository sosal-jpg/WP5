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
t_col=0.003
t_sand=0.003
F_z=1e3
F_x,F_y=1e3

def mass(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D
    return mat_att.rho(((a+b+t)*w-4*pi*D**2/4)*t)+4*mat_fast.rho*(t_col + t + 0.8 * D * 2 * 1.5 ** 2)*(pi*D**2/4)
def sandwich_side_failure(x):
    a_,b_,t,w_,D=x
    return ((F_xy/(D*t))**2+3*)**(1/2)

def lug_hole_sand(x):
    a_,b_,t,w_,D=x
    s_plane=(F_y**2+F_x**2)**(1/2)/(t*D)
    s_shear_z=F_z/(pi*D*t)
    return (s_plane**2+3*s_shear_z)**(1/2)

def sand_hole_sand(x):
    a_,b_,t,w_,D=x
    s_plane=(F_y**2+F_x**2)**(1/2)/(t_sand*D)
    s_shear_z=F_z/(pi*D*t_sand)
    return (s_plane**2+3*s_shear_z)**(1/2)

def bolt_hole_sand(x):
    a_,b_,t,w_,D=x
    A=pi*D**2/4
    s_shear=(F_y**2+F_x**2)**(1/2)/A
    s_axial=F_z/A
    return (s_axial**2+3*s_shear)**(1/2)

def lug_hole_col(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D
    F_yz=F_z*(a-e+t/2)/(b-e-t)
    F_yx=F_x*(a-e+t/2)/(w-2*e)
    F_zx=F_x*(b-e+t/2)/(w-2*e)
    s_plane=((F_z+F_zx)**2+F_x**2)**(1/2)/(t*D)
    s_shear_z=(F_y+F_yz+F_yx)/(pi*D*t)
    return (s_plane**2+3*s_shear_z)**(1/2)

def col_hole_col(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D
    F_yz=F_z*(a-e+t/2)/(b-e-t)
    F_yx=F_x*(a-e+t/2)/(w-2*e)
    F_zx=F_x*(b-e+t/2)/(w-2*e)
    s_plane=((F_z+F_zx)**2+F_x**2)**(1/2)/(t_col*D)
    s_shear_z=(F_y+F_yz+F_yx)/(pi*D*t_col)
    return (s_plane**2+3*s_shear_z)**(1/2)

def bolt_hole_col(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D
    A=pi*D**2/4
    F_yz=F_z*(a-e+t/2)/(b-e-t)
    F_yx=F_x*(a-e+t/2)/(w-2*e)
    F_zx=F_x*(b-e+t/2)/(w-2*e)
    s_shear=((F_z+F_zx)**2+F_x**2)**(1/2)/A
    s_axial=(F_y+F_yz+F_yx)/A
    return (s_axial**2+3*s_shear)**(1/2)

def bending(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D

