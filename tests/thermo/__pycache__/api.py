
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
        with open(coord_file, 'r') as f:
            line = f.readline()
            line = line.split()
            data_N = int(line[0])
            cell = np.array([float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6]),])
            if cell == None:
                pbc = False
            else:
                pbc = True
            line = f.readline()
            data_atoms = np.array([])
            data_xyz = np.zeros((data_N, 3))
            for i in range(0, data_N ):
                line = f.readline()
                line = line.split()
                data_atoms = np.append(str(data_atoms), line[0])
                data_xyz[i, 0] = float(line[1])
                data_xyz[i, 1] = float(line[2])
                data_xyz[i, 2] = float(line[3])


        return data_N, data_atoms, data_xyz, cell, pbc
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
            line = f.readline()
            line = line.split()
            data_N = int(line[0])
            system = str(line[1])
            if system == "S":
                pbc = True
            else:
                pbc = False 
            line = f.readline()
            atomic_species = str(line.split())
            data_xyz = np.zeros((data_N, 3))
            data_atoms = np.array([])
            for i in range(0, data_N ):
                line = f.readline()
                line = line.split()
                index = int(line[0])
                symbol_number = int(line[1])
                data_xyz[i, 0] = float(line[2])
                data_xyz[i, 1] = float(line[3])
                data_xyz[i, 2] = float(line[4])
                data_atoms = np.append(str(data_atoms), atomic_species[symbol_number-1])
            if pbc:
                line = f.readline()
                line = line.split()
                cell_vector = np.array([[float(line[0]), float(line[1]), float(line[2])], [float(line[3]), float(line[4]), float(line[5])], [float(line[6]), float(line[7]), float(line[8])]])
            else:
                cell_vector = None

        return data_N, data_atoms, data_xyz, cell_vector, pbc
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
        line = f.readline()
        line = line.split()
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
        elif coord_file.endswith(".gen"):
            data_N, data_atoms, data_xyz, cell, pbc = read_gen(coord_file)
        else:
            SystemError("The input file is not supported.")
        
    else:
        SystemError("The engine is not supported.")

    return [data_N,data_atoms,data_xyz,cell]

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
    system_info = System(atoms=atom_list,electronic_energy=energy,cell=cell,vibrational_frequencies=vibrational_frequencies)
    print("All values of system_info:")
    print("############################################### ")
    print("Atoms:", system_info._atoms)
    print("Electronic energy: ", system_info._electronic_energy)
    print("Cell: ", system_info._cell)
    print("Vibrational frequencies: ", system_info._vibrational_frequencies)
    print("Number of atoms: ", system_info.number_of_atoms)
    print("Dimension: ", system_info.dim)
    print("Degree of freedom: ", system_info.dof)
    print("Charge: ", system_info.charge)
    print("Mass: ", system_info.mass)
    print("Spin: ", system_info.spin)
    print("Rotational symmetry number: ", system_info.rotational_symmetry_number)
    print("Spacegroup number: ", system_info.spacegroup_number)
    print("Spacegroup: ", system_info.spacegroup)
    print("Periodicity: ", system_info.periodicity)
    print("PBC: ", system_info.pbc)
    print("Solvation: ", system_info.solvation)
    print("Solvent: ", system_info.solvent)
    print("Center of mass: ", system_info.center_of_mass)
    print("############################################### ")
    print("System object initialized")
    print("-----------------------------------------------")

    # Initialize thermo
    print("Initializing Thermo object")
    thermo_setup = Thermo(system=system_info,temperature=temperature,pressure=pressure)
    print("Thermo object initialized")
    print("-----------------------------------------------")
    # Calculate thermo
    print("Calculating thermo")
    thermo_setup.run()
    print("Thermo calculated")
    print("-----------------------------------------------")

  

    






    



