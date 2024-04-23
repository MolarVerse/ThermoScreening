
from .inputFileReader import InputFileReader
from .system import System
from .thermo import Thermo
from .atoms import Atom
import numpy as np

def read_xyz(coord_file: str):
        """
        Reads xyz-file and returns number of atoms, chemical symbols of atoms, coordinates of atoms, and cell parameters.

        Parameter:
        ----------
        coord_file : str
            The name of the coordinate file.
        
        Returns:
        --------
        List
            number of atoms, chemical symbols of atoms, coordinates of atoms, cell parameters
        """
        
  
        cell = None
        pbc = False
        data_N = None
      
        with open(coord_file, 'r') as f:
            
            line = f.readline()
 
            line = line.strip()
           
            # split line and ignore spaces
            line = line.split()
            data_N = int(line[0])

           
            # if len(line) > 1 and line[1:] is containing only numbers
            if len(line) > 1 and all([x.replace('.','',1).isdigit() for x in line[1:]]):
                
                cell = np.array([float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6])])
                pbc = True

            line = f.readline().strip()
     
            i = 0
            data_atoms = np.zeros(data_N, dtype=str)
            data_xyz = np.zeros((data_N, 3), dtype=float)
            while True:
                line = f.readline().strip()
                line = line.split()
                if len(line) == 0:
                    break
                data_atoms[i] = str(line[0])
                data_xyz[i, 0] = float(line[1])
                data_xyz[i, 1] = float(line[2])
                data_xyz[i, 2] = float(line[3])
                i += 1
             



        return [data_N, data_atoms, data_xyz, cell, pbc]
def read_gen(coord_file: str):
        """
        Reads gen-file and returns number of atoms, chemical symbols of atoms, coordinates of atoms, and cell vector.

        Parameter:
        ----------
        coord_file : str
            The name of the coordinate file.
        
        Returns:
        --------
        List
            number of atoms, chemical symbols of atoms, coordinates of atoms, cell vector
        """
        with open(coord_file, 'r') as f:
            line = f.readline().strip()
            line = line.split()
            data_N = int(line[0])
            system = str(line[1])
            if system == "S":
                pbc = True
            else:
                pbc = False 
            line =f.readline().strip()
            atomic_species = str(line.split())
            data_xyz = np.zeros((data_N, 3))
            data_atoms = np.array([])
            i = 0
            
            data_atoms = np.zeros(data_N, dtype=str)
            data_xyz = np.zeros((data_N, 3), dtype=float)
            while True:
                line = f.readline().strip()
                line = line.split()
                if i < data_N:
                    break
                index = int(line[0])
                symbol_number = int(line[1])
                data_xyz[i, 0] = float(line[2])
                data_xyz[i, 1] = float(line[3])
                data_xyz[i, 2] = float(line[4])
                data_atoms = np.append(str(data_atoms), atomic_species[symbol_number-1])
                i += 1
            if pbc:
                line = f.readline()
                line = f.readline()
                line = line.split()
                line2 = f.readline()
                line2 = line2.split()
                line3 = f.readline()
                line3 = line3.split()
                cell_vector = np.array([[float(line[0]), float(line[1]), float(line[2])], [float(line2[0]), float(line2[1]), float(line2[2])], [float(line3[0]), float(line3[1]), float(line3[2])]])
            else: cell_vector = None

        return [data_N, data_atoms, data_xyz, cell_vector, pbc]
def read_vib_file(vibrational_file: str):
    """
    Reads the vibrational file and returns the vibrational frequencies as a numpy array.

    Parameter:
    ----------
    vibrational_file : str
        The name of the vibrational file.
    
    Returns:
    --------

    np.ndarray
        The vibrational frequencies as a numpy array.
    """
    vibrational_frequencies = np.array([])
    with open(vibrational_file, 'r') as f:
        while True:
            line = f.readline().strip()
            line = line.split()
            if len(line) == 0:
                break
            vibrational_frequencies = np.append(vibrational_frequencies, float(line[1]))
         

    return vibrational_frequencies

def read_coord(coord_file: str,engine: str):
    """
    Reads the coordinate file and returns the coordinates as a numpy array.

    Parameters
    ----------
    coord_file : str
        The name of the coordinate file.
    engine : str
        The name of the engine.
    
    Returns
    -------
    List
        number of atoms, chemical symbols of atoms, coordinates of atoms
    """
    print("coord_file: ", coord_file)
    print("engine: ", engine)
    data_N, data_atoms, data_xyz, cell = None, None, None, None
    if engine == "dftb+":
        if coord_file.endswith(".xyz"):
            data_N, data_atoms, data_xyz, cell, pbc = read_xyz(coord_file)
            return data_N,data_atoms,data_xyz,cell,pbc
        elif coord_file.endswith(".gen"):
            # data_N, data_atoms, data_xyz, cell_vectors, pbc = read_gen(coord_file)
            # return data_N,data_atoms,data_xyz,cell_vectors,pbc
            raise SystemError("The gen file is not tested yet.")
        else:
            SystemError("The input file is not supported.")
        
    else:
        SystemError("The engine is not supported.")

    

def read_vibrational(vibrational_file: str,engine: str):
    """
    Reads the vibrational file and returns the vibrational frequencies as a numpy array.
    
    Parameters
    ----------
    vibrational_file : str
        The name of the vibrational file.
    engine : str

    Returns
    -------
    np.ndarray
        The vibrational frequencies as a numpy array.
    """
    vibrational_frequencies = None
    if engine == "dftb+":
        vibrational_frequencies = read_vib_file(vibrational_file)
    else:
        # throw error
        SystemError("The engine is not supported.")
   
    return vibrational_frequencies
def unit_length(engine : str):
    """
    Returns the unit of length.

    Parameter:
    ----------
    engine : str
        The name of the engine.
    
    Returns:
    --------
    str
        The unit of length.
    """
    if engine == "dftb+":
        return "Angstrom"
    else:
        SystemError("The engine is not supported.")
def unit_energy(engine : str):
    """
    Returns the unit of energy.

    Parameter:
    ----------
    engine : str
        The name of the engine.
    
    Returns:
    --------
    str
        The unit of energy.
    """
    if engine == "dftb+":
        return "Hartree"
    else:
        SystemError("The engine is not supported.")
def unit_mass(engine : str):
    """
    Returns the unit of mass.

    Parameter:
    ----------
    engine : str
        The name of the engine.
    
    Returns:
    --------
    str
        The unit of mass.
    """
    if engine == "dftb+":
        return "amu"
    else:
        SystemError("The engine is not supported.")     
def unit_frequency(engine : str):
    """
    Returns the unit of frequency.

    Parameter:
    ----------
    engine : str
        The name of the engine.
    
    Returns:
    --------
    str
        The unit of frequency.
    """
    if engine == "dftb+":
        return "cm^-1"
    else:
        SystemError("The engine is not supported.")


def execute(input_file: str):
    # Read input file
    print("-----------------------------------------------")
    print("Reading input file")
    input_file_reader = InputFileReader(input_file)
    input_dictionary = input_file_reader._dictionary

   
    print("Input file read")
    print("-----------------------------------------------")
    # Read coordinate file
    print("Reading coordinate file")
    coordinate_file = input_dictionary["coord_file"]
    engine = input_dictionary["engine"]
    atom_number,atom_symbol,xyz,cell,pbc = read_coord(coordinate_file,engine)
    # print("xyz ", xyz)
    # print("atom_symbol ", atom_symbol)
    print("Coordinate file read")
    print("-----------------------------------------------")

    # Read vibrational file
    print("Reading vibrational file")
    vibrational_file = input_dictionary["vibrational_file"]
    vibrational_frequencies = read_vibrational(vibrational_file,engine)
    print("Vibrational file read")
    print("-----------------------------------------------")

    # Read temperature
    print("Reading temperature")
    temperature = float(input_dictionary["temperature"])
    print("Temperature read: T = ",temperature, " K")
    print("-----------------------------------------------")

    # Read pressure
    print("Reading pressure")
    pressure = float(input_dictionary["pressure"])
    print("Pressure read: p = ",pressure, " Pa")
    print("-----------------------------------------------")

    # Read energy
    print("Reading energy")
    energy = float(input_dictionary["energy"])
    print("Energy read: E = ",energy, " H")
    print("-----------------------------------------------")

    # Initialize atoms
    atom_list = []
    print("Initializing Atom objects")
    for i,symbol in enumerate(atom_symbol):
     
        atom_list.append(Atom(symbol=symbol,position=xyz[i,:]))
    
    print("Atom objects initialized")
    print("-----------------------------------------------")

   
 

    # Initialize system
    print("Initializing System object")
    system_info = System(atoms=atom_list,electronic_energy=energy,cell=cell,vibrational_frequencies=vibrational_frequencies,pbc=pbc)
    # check vibrational_frequencies
    print("#################################################################### ")
    if system_info._check_frequency_length == False:
        ValueError("The number of vibrational frequencies does not match with the degree of freedom.")
    print("##################### All values of system_info: ######################")
    print("#####################################################################")
    print("Atoms:", system_info.atom_names())
    print("Coordinates: ", system_info.coord(), "in ",unit_length(engine)) 
    print("Electronic energy: ", system_info._electronic_energy, " in ", unit_energy(engine))
    print("Cell: ", system_info._cell, " in ", unit_length(engine))
    print("Vibrational frequencies: ", system_info._vibrational_frequencies, " in ", unit_frequency(engine))
    print("Real vibrational frequencies: ", system_info._real_vibrational_frequencies, " in ", unit_frequency(engine))

    print("Number of atoms: ", system_info.number_of_atoms)
    print("Dimension: ", system_info.dim)
    print("Degree of freedom: ", system_info.dof)
    print("Charge: ", system_info.charge)
    print("Mass: ", system_info.mass, " in ", unit_mass(engine))
    print("Spin: ", system_info.spin)
    print("Rotational symmetry number: ", system_info.rotational_symmetry_number())
    print("Rotational Group: ", system_info.rotational_group)
    # print("Spacegroup number: ", system_info.spacegroup_number())
    # print("Spacegroup: ", system_info.spacegroup())
    print("Periodicity: ", system_info.periodicity)
    print("PBC: ", system_info.pbc)
    print("Solvation: ", system_info.solvation)
    print("Solvent: ", system_info.solvent)
    print("Center of mass: ", system_info.center_of_mass , " in ", unit_length(engine))
    print("############################################### ")
    print("System object initialized")
    print("-----------------------------------------------")

    # Initialize thermo
    print("Initializing Thermo object")
    thermo_setup = Thermo(system=system_info,temperature=temperature,pressure=pressure,engine=engine)
    print("Thermo object initialized")
    print("-----------------------------------------------")
    # Calculate thermo
    print("Calculating thermo")
    thermo_setup.run()
    print("Thermo calculated")
    print("-----------------------------------------------")

  

    






    



