# -*- coding: utf-8 -*-
"""
@author: Stefanie Kroell
"""
### Import used libaries ###

import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt

### Define physcical constants ###

h_quer = np.double(const.hbar)
pi = np.double(const.pi)
c = np.double(const.c)
kB = np.double(const.k)
R = np.double(const.gas_constant)
h = np.double(const.h)
u = np.double(const.u)
kcal = const.calorie/1000
angstrom = np.double(const.angstrom)
eV = np.double(const.electron_volt)
H = const.physical_constants["Hartree energy"][0]
amu = const.physical_constants["atomic unit of mass"][0]
bohr_radius = 5.29177210903*10**(-1)
Na = const.N_A
cal = const.calorie



### Define temperature and pressure of the system ###

T = 298.15 # in unit K
P = 101325 # in unit Pa

### User input ###

#sigma_r = 0 

# symmetry number for rotation by manual user input

sigma_r = int(input("Enter the symmetry number for rotation sigma_r: "))

#s = 0  #2s+1

# Manual user input for electronic spin multiplicity

s = int(input("Enter the electronic spin multiplicity of the molecule s: "))

#dof = 6

# Manual user input of degree of freedom

dof = int(input("If it is a linear Molecule enter 5, if it is not enter 6: "))


#E0 = -155.046208186

### Import the calculated xyz-datafile for geomety and atom type ###

path = os.getcwd()
path_molecule = path+"/geo_opt.xyz"
# number of atoms
data_N = genfromtxt(path_molecule, delimiter='   ', \
usecols=0,max_rows=1) 
data_atoms = genfromtxt(path_molecule, delimiter='   ', \
skip_header=2,usecols=0,dtype=str) # string
data_xyz = genfromtxt(path_molecule, delimiter='     ', \
skip_header=2,usecols=(1,2,3)) # in unit Angstrom

### Import the calculated electronic energy ###

pathE0 = path+"/electronic_energy.txt"
E0 = genfromtxt(pathE0,delimiter=' ', usecols=0, max_rows=1)
print("Total Energy: {} in H".format(E0)) # in unit Hartree
print("\n")

### Begin of thermochemical calculation ###

def start():
    

    # column of the atom name
    atoms_name = data_atoms[:]

    #the coordinates converted into metre

    # angstrom = 10**(-10)
    # assign the columns of the data file to coordinates
    x = data_xyz[:,0] # in angstrom 
    y = data_xyz[:,1] # in angstrom 
    z = data_xyz[:,2] # in angstrom

    # check the import
    print("Atomsorten: {}".format(atoms_name))
    print("\n")
    print("x-Koordinaten in Angstrom: {}".format(x))
    print("\n")
    print("y-Koordinaten in Angstrom: {}".format(y))
    print("\n")
    print("z-Koordinaten in Angstrom: {}".format(z))
    print("\n")

    # check the number of atoms

    N = len(atoms_name)

    # check the number of atoms

    print("Anzahl der Atome: {}".format(data_N))
    print("Anzahl der Zeilen: {}".format(N))
    print("\n")

    ### Import list of atom mass in Molar mass unit ###
    path_mass = path+"/mass.txt"
    m_atom_name = genfromtxt(path_mass, delimiter='	', \
    skip_header=2,usecols=0,dtype=str)
    m_u = genfromtxt(path_mass, delimiter='	',skip_header=2,usecols=1)

    # convert into kg
    # m = m_u*u #in unit kg

    # not converted
    m = m_u

    # create an numpy array of atom masses
    mass_atom = np.array([m_atom_name,m])

    # check
    print("Masse der Atome: {} in u".format(mass_atom))
    print("\n")



    # assign the sort of atom to the mass

    def mass_assignment(atoms_name,N,mass_atom):
        mass_pos = np.empty(0) # create an empty list mass_pos
        # length of the data with the mass of the atomsorts
        N_atomsorts = len(mass_atom[0][:]) 

        for i in range(N):
            for j in range(N_atomsorts):
                #check if the atom name of the xyz file is the same name of 
                # the atom name of the data file with the mass
                if atoms_name[i] == mass_atom[0][j]: 
                    mass_pos = np.append(mass_pos,mass_atom[1][j]) 
                    # if it is truth then on that position the assosciated
                    # mass is saved into the list mass_pos
        # return the converted sting to the function
        return(np.double(mass_pos)) 
       

    # the associated mass of the xyz position
    m_pos =  mass_assignment(atoms_name,N,mass_atom) 

    # check
    print("Zudordnung der Massen (in u) zu den Atomen:{}".format(m_pos))
    print("\n")
    M = np.sum(m_pos)
    print("Summe Masse (in u):{}".format(M))
    print("\n")

    ### Compute mass center position ###
    def mass_center_position(m_pos,x,y,z,N):
        # Sum over atoms
        def sum_mr(m,r,N):
            mr = np.empty(0)
            for i in range(N):
                mr = np.append(mr,m[i]*r[i])
            #print("Summe Masse und Ort:{}".format(mr))
            # Total mass
            sum_mr = np.sum(mr)
            # print(sum_mr)
            return(sum_mr)

        # for each direction
        sum_mx = sum_mr(m_pos,x,N)
        sum_my = sum_mr(m_pos,y,N)
        sum_mz = sum_mr(m_pos,z,N)

        # print(sum_mx)
        # print(sum_my)
        # print(sum_mz)

        Rx = sum_mx/M
        Ry= sum_my/M
        Rz = sum_mz/M

        R = [Rx,Ry,Rz]
        # print(R)
        # return the mass point 
        return(R)

    # call t he mass_center_position function and save return value
    #rsp = mass_center_position(m_pos,x,y,z,N)
    #rsp =
    #print(rsp)

    ### relocate the coordinates relative to the center of mass R ###

    def relocate(x,y,z,rsp):
        rx = np.subtract(x,rsp[0])
        ry = np.subtract(y,rsp[1])
        rz = np.subtract(z,rsp[2])

        print(rx,"\n", ry,"\n",rz)
        r = [rx,ry,rz]
        return(r)

    # transformed coordinates
    #relR = relocate(x,y,z,rsp)
    relR = [x,y,z]

    #print("Koordinaten im Schwerpunkt:{}".format(relR))

    ### calculate the inertia tensor ###
    def inertia_tensor(relR,m_pos,N):
        x = relR[0]
        y = relR[1]
        z = relR[2]
        ### diagonal intertia tensor elements ###
        def Iii_function(ri,rj,m_pos,N):
            i_ii = np.empty(0)
            for i in range(N):
                mii = m_pos[i]*(ri[i]**2+rj[i]**2)
                i_ii = np.append(i_ii, mii)
            #print("i_ii:{}".format(i_ii))
            Iii = np.sum(i_ii)
            #print("Summe Iii:{}".format(Iii))
            return(Iii)
        ### non-diagonal intertia tensor elements ###
        def Iij_function(ri, rj, m_pos, N):
            i_ij = np.empty(0)
            for i in range(N):
                mij = m_pos[i] * (ri[i]*rj[i])
                i_ij = np.append(i_ij, mij)
            #print("i_ij:{}".format(i_ij))
            Iij = -np.sum(i_ij)
            #print("Summe Iij:{}".format(Iij))
            return(Iij)
        ### all tensor elements ###
        Ixx = Iii_function(y, z ,m_pos,N)
        Iyy = Iii_function(x, z, m_pos, N)
        Izz = Iii_function(x, y, m_pos, N)
        Ixy = Iij_function(x, y, m_pos, N)
        Iyz = Iij_function(y, z, m_pos, N)
        Ixz = Iij_function(x, z, m_pos, N)

        #print(Ixx)

        ### Store tensor in array ###

        I = np.array([[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]])

        # check
        print("Massenträgheitstensor in u °A²: \n{}".format(I))
        print("\n")
        # return inertia tensor array
        return(I)
    # Inertia tensor initialization
    I = inertia_tensor(relR,m_pos,N)

    ### calculate the eigenvecotrs and eigenvalues ###

    eigenvalues = np.linalg.eig(I) # in unit u*A^2
    # eigenvalue
    evl = np.abs(eigenvalues[0])
    # eigenvector
    evv = eigenvalues[1]

    print("Eigenwerte und Eigenvektoren in u °A²: {}".format(eigenvalues))
    print("\n")
    print("Eigenwert berechnet in u °A²:{}".format(evl)) 
    print("Eigenvektor in u °A²: {}".format(evv))
    print("\n")


    ##############################################
    ### Computional of rotational contribution ###
    ##############################################

    # Transform eigenvalue in SI units
    evl_SI = evl*u*(10**(-20)) # in kg*m^2

    # Rotational temperatures in x-,y-,z-direction
    rot_temp = (((h**2)/(8*pi**2*kB))/evl_SI) # in  K
    # Rotational constants in x-,y-,z-direction
    rot_const = ((h_quer**2/2)/evl_SI)/h*10**(-9) # in GHz

    # Total rotational temperature
    rot_temp_xyz = rot_temp[0]*rot_temp[1]*rot_temp[2] # in K
    # Check
    print("Rot_temp_xyz: {}".format(rot_temp_xyz))
    
    # partition function contribution
    q_r = (pi**(1/2)/sigma_r)*(T**(3/2)/(np.power(rot_temp_xyz,1/2)))

    # entropy contribution
    S_r = R*(np.log(q_r)+3/2)/cal # in cal/mol K

    # energy contribution
    E_r = ((3/2)*R*T)/H/Na # in Hartree/Particle
    E_r_cal = E_r*H*Na/cal # in cal

    # heat capacity contribution
    C_r = (3/2)*R/cal # in cal/ mol K

    # Check
    print("Rotationsanteil: Rot_temp in K:{}, \
     Rot_const in GHz:{}".format(rot_temp,rot_const))
    print("\n")	
    print("q_r: {}, E_r: {} in H/Particle, S_r: {} in cal/(mol K),\
     C_r: {} in cal/(mol K)".format(q_r,E_r,S_r,C_r))
    print("\n")
    
    ### vibrational contribution (BOT) ###

    # Input computed frequencies ###
    path_nu = path + "/frequency.txt"
    nu_K = genfromtxt(path_nu, delimiter='  ',usecols=2,skip_header=dof)

    # Check
    print("Frequenzen in cm⁻1:{}".format(nu_K))
    print("\n")

    # Vibrational temperature in x-,y-,z-direction
    vib_temp_K = h*(nu_K*c*10**(2))/kB # in K
    

    
    # total vibrational temperature
    s_v = np.sum(vib_temp_K/2) # in K

    # Check
    print("Summe vib_temp_K:{}".format(s_v))
    
    # partition function contribution element
    q_v_K = np.double(np.exp(-vib_temp_K/(2*T))/(1-np.exp(-vib_temp_K/T)))

    # product of q_v_K
    q_v = np.prod(q_v_K)

    # entropy contribution element
    S_v_sum = np.subtract(np.divide((vib_temp_K/T),(np.exp(vib_temp_K/T)-1)),\
    np.log(1-np.exp(-vib_temp_K/T))) # in cal/mol K
    
    # Check
    #print("S_v_sum:{}, len:{}".format(S_v_sum,len(S_v_sum)))

    # entropy contribution
    S_v = R*np.sum(S_v_sum)/cal

    # zero point energy correction 
    E_corr= R*(np.sum(vib_temp_K/2))/H/Na # in H

    # energy contribution
    # in Hartree/Mol
    E_v = R*np.sum(np.multiply(vib_temp_K,(1/2+1/(np.exp(vib_temp_K/T)-1))))/H/Na 
    E_v_cal = E_v*H*Na/cal # in cal

    # C_v = (R*np.sum(np.multiply(np.exp(vib_temp_K/T), \
    # (np.power(np.divide((vib_temp_K/T),(np.exp(-vib_temp_K/T)-1)), 2)))))
    #vib_temp_K_C = h*(nu_K*c)/kB # in K
    
    # heat capacity contriubtion
    C_v = R*np.sum(np.exp(-vib_temp_K/T)*((vib_temp_K/T)/ \
    (np.exp(-vib_temp_K/T)-1))**2)/cal # in in cal/ mol K
    EZP = R*(np.sum(vib_temp_K/2))/cal #in cal
    
    # Check
    print("vib_temp: {} in K,len:{}".format(vib_temp_K,len(vib_temp_K)))
    print("\n")
    print("q_v: {}, E_v: {} in H/Particle, S_v: {} in cal/(mol K), \
    C_v: {} in cal/(mol K)".format(q_v,E_v,S_v, C_v))
    print("\n")
    print("Zero Point correction: {} in H/Particle".format(E_corr))
    print("\n")
    print("Zero Point vibrational energy: {} in cal/mol".format(EZP))
    print("\n")
    
    ### translation contribution ###
    m_M = np.sum(m_pos)*u #Molecule mass in kg
    
    #print("m_M:{}".format(m_M))
    
    # partition function contribution
    q_t = np.power((2*pi*m_M*kB*T)/(h**2),(3/2))*(kB*T/P) 

    # entropy contribution
    S_t = R*(np.log(q_t)+5/2)/cal # in cal / mol K

    # energy contribution
    E_t = (3/2*R*T)/H/Na # in H/Particle
    E_t_cal = E_t*H*Na/cal # in cal

    # heat capacity contribution 
    C_t = 3/2*R/cal # in cal/ mol K

    # Check
    print("q_t: {}, E_t: {} in H/Particle, S_t: {} in cal/(mol K),  \
     C_t: {} in cal/(mol K)".format(q_t,E_t,S_t,C_t))
    print("\n")


    ### electronic contribution ###

    # spin multiplicity
    omega_0 = 2*s+1 # s=> spin: user input

    # partition function contribution
    q_e = omega_0

    # energy contribution
    E_e = 0/H/Na

    # entropy contribution
    S_e = R*(np.log(q_e))/cal # in cal
    
    # heat capacity contribution
    C_e = 0

    # Check
    print("q_e: {}, E_e: {} in H/Particle, S_e: {} in cal/(mol K), \
     C_e: {} in cal/(mol K)".format(q_e,E_e,S_e, C_e))

    ### total internal thermal energy ###

    E_tot = E_t+E_r+E_v+E_e #thermal correction to energy
    E_tot_cal = E_tot*H*Na/cal # in cal

    # Check
    print("E_corr: {} in H/particle, \
     E_corr_cal:{} in cal/mol".format(E_tot,E_tot_cal))
    print("\n")

    ### total entropy ###

    S_tot = (S_t+S_r+S_v+S_e) # in cal / mol K

    # Check 
    print("S_corr: {} in cal/(mol K)".format(S_tot))
    print("\n")

    ### Total enthalpy ###

    H_corr = E_tot+(R*T)/H/Na # in Hartree/Particle

    # Check 
    print("H_corr: {} in H/Particle".format(H_corr))
    print("\n")


    ### Gibss Free energy ###

    G_corr = H_corr-(S_tot*cal/Na/H)*T # in Hartree/Particle

    # Check 
    print("G_corr: {} in H/Particle".format(G_corr))
    print("\n")


    ### heat capacity ###

    C_tot = (C_t+C_r+C_v+C_e) #in cal/mol K

    # Check
    print("C_corr: {} in cal/(mol K)".format(C_tot))
    print("\n")

    ### Summmarize ###

    EeZP = (E0+E_corr)

    print("Sum of electronic and zero-point Energies: \
    {} in H/Particle".format(EeZP))

    EeEtot = (E0+E_tot)

    print("Sum of electronic and thermal Energies: \
    {} in H/Particle".format(EeEtot))
    print("\n")
    
    EeHcorr = (E0+H_corr)

    print("Sum of electronic and thermal Enthalpies: \
    {} in H/Particle".format(EeHcorr))
    print("\n")
     
    EeGcorr = (E0+G_corr)

    print("Sum of electronic and thermal Free Energies: \
    {} in H/Particle".format(EeGcorr))
    print("\n")
    
	
start()
