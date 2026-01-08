import numpy as np
from numpy import sin, cos,pi
import math
from scipy.optimize import minimize, NonlinearConstraint

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

mat_att=Material(Name='7075-T6 Aluminium', E=73.1e9, y_stress=503e6, alpha=23.6e-6, rho=2810)
mat_fast=Material(Name='7075-T6 Aluminium', E=73.1e9, y_stress=503e6, alpha=23.6e-6, rho=2810)
mat_sand=Material(Name='7075-T6 Aluminium', E=73.1e9, y_stress=503e6, alpha=23.6e-6, rho=2810)
mat_col=Material(Name='8552 Carbon-Epoxy (Hexcel)', E=190e9, y_stress=3310e6, alpha=6e-6, rho=1300)
t_col=0.0004
t_sand=0.0002

g=9.8
sandwich=[[6.27+2,0],[37+2,0.6],[7.77+2,1.2],[2,1.8]]
axial_a=6*g
transterse_a=2*g
n=4
F_y=F_x=transterse_a*sandwich[1][0]
F_z=axial_a*sandwich[1][0]



def mass(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D
    return mat_att.rho*(((a+b+t)*w-4*pi*D**2/4)*t)+4*mat_fast.rho*(0.0016 + t + 0.8 * D * 2 * 1.5 ** 2)*(pi*D**2/4)

def lug_hole_sand(x):
    a_,b_,t,w_,D=x
    s_plane=(F_y**2+F_x**2)**(1/2)/(t*D)
    s_shear_z=F_z/(pi*D*t)
    return (s_plane**2+3*s_shear_z)**(1/2)/mat_att.y_stress/2

# def sand_hole_sand(x):
#     a_,b_,t,w_,D=x
#     s_plane=(F_y**2+F_x**2)**(1/2)/(t_sand*D)
#     s_shear_z=F_z/(pi*D*t_sand)
#     return (s_plane**2+3*s_shear_z)**(1/2)/mat_sand.y_stress

def sand_hole_sand(x):
    a_,b_,t,w_,D=x
    s_plane=0
    s_shear_z=F_z/(pi*D*t_sand)
    return (s_plane**2+3*s_shear_z)**(1/2)/mat_sand.y_stress/2

def bolt_hole_sand(x):
    a_,b_,t,w_,D=x
    A=pi*D**2/4
    s_shear=(F_y**2+F_x**2)**(1/2)/A
    s_axial=F_z/A
    return (s_axial**2+3*s_shear)**(1/2)/mat_fast.y_stress/2

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
    return (s_plane**2+3*s_shear_z)**(1/2)/mat_att.y_stress/2

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
    return (s_plane**2+3*s_shear_z)**(1/2)/mat_col.y_stress/2

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
    return (s_axial**2+3*s_shear)**(1/2)/mat_fast.y_stress/2

def bending(x):
    a_,b_,t,w_,D=x
    e=1.5*D
    a=a_+2*e
    b=b_+2*e
    w=w_+2*e+3*D
    A=w*t
    I_y=t*w**3/12
    R_c=1.5*t
    R_n=t/math.log(2)
    e=R_c-R_n
    M_x=F_z*(a-e+t/2)+F_y*(1.5*t)
    s_z=M_x/(A*e)*(R_n/t-1)+F_x*(1.5*t)*w/2/I_y+F_z/A
    return s_z/mat_att.y_stress

constraints=[
    NonlinearConstraint(lug_hole_sand, lb=0,ub=1),
    NonlinearConstraint(sand_hole_sand, lb=0,ub=1),
    NonlinearConstraint(bolt_hole_sand, lb=0,ub=1),
    NonlinearConstraint(lug_hole_col, lb=0,ub=1),
    NonlinearConstraint(bolt_hole_col, lb=0,ub=1),
    NonlinearConstraint(bending, lb=0,ub=1),
]
options = {"maxiter": 1000, }
x0=[0,0,0.001,0,0.001]
bounds=[
    (0,np.inf),
    (0,np.inf),
    (0.001,np.inf),
    (0,np.inf),
    (0.001,np.inf),
]

res = minimize(mass, x0=x0, constraints=constraints, bounds=bounds, method='SLSQP', options=options)
print(res)
for c in constraints:
    print(c.fun(res.x))
print(col_hole_col(res.x))

