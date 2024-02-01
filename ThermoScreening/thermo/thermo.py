import numpy as np
from .physicalConstants import PhysicalConstants
from .system import System


class Thermo:
    """
    A class for calculating the thermochemistral properties e.g. Gibbs free energy, enthalpy, entropy, heat capacity, etc. of a system.
    The class needs as input the temperature, pressure, system information (coordinates, degree of freedom, electronic spin, symmetry number, electronic energy and vibrational frequencies).
    """

    def __init__(self, temperature: float, pressure: float,system: System,engine: str):
        """
        Parameters
        ----------
        temperature : float
            The temperature of the system in Kelvin.
        pressure : float
            The pressure of the system in bar.
        system : System
            The system information of the system.
        engine : str
            The engine used for the calculation to compute the thermochemical properties with correct units.
        """
        if pressure == None:
           ValueError("The pressure is not given.")
        if temperature == None:
            ValueError("The temperature is not given.")   
        if system == None:
            ValueError("The system is not given.")
        if engine == None:
            ValueError("The engine is not given.")
        self._temperature = temperature
        self._pressure = pressure
        self._system = system
        self._engine = engine
        if self._engine != "DFTB+":
           ValueError("The engine is not supported.")
        if self._temperature < 0:
           ValueError("The temperature is negative.")
        if self._pressure < 0:
           ValueError("The pressure is negative.")
        self._coord = self._system.coord()
        self._atomic_masses = self._system.atomic_masses()


    def run(self):
        """
        Calculates the thermochemical properties of the system.
        """
        #self._transform_units()
        self._rotational_contribution()
        self._vibrational_contribution()
        self._electronic_contribution()
        self._translational_contribution()
        
        self._compute_thermochemical_properties()
        self._summary()

    def _transform_units(self):
        if self._engine == "DFTB+":
            # the originial units of DFTB+ are Hartree, Angstrom, amu
            # the units are transformed to J, m, kg
            self._coord = self._system.coord()[:]*PhysicalConstants["A"]
            self._atomic_masses = self._system.atomic_masses()[:]*PhysicalConstants["u"]
            self._system.real_vibrational_frequencies[:] = self._system.real_vibrational_frequencies[:]*PhysicalConstants["HztoGHz"]
            self._system.electronic_energy = self._system.electronic_energy*PhysicalConstants["H"]*PhysicalConstants["Na"]
            self._system.real_vibrational_frequencies[:] = self._system.real_vibrational_frequencies[:]*PhysicalConstants["HztoGHz"]
      
         
        else:
            ValueError("The engine is not supported.")

    def _relocate_to_cm(self):
        """
        Relocates the system to the center of mass.
        """
 
        coord = self._system.coord()
        coord -= self._system.center_of_mass
        # print("########################################")
        # print("Relocate to center of mass: ")
        # print(coord)
    def _compute_inertia_tensor(self):
        """
        Computes the inertia tensor of the system.
        """
        coord = self._system.coord()
        x = coord[:,0]
        y = coord[:,1]
        z = coord[:,2]
  
        m = self._system.atomic_masses()
        
       
        
        
        self._inertia_tensor = np.zeros((3,3))

        self._inertia_tensor[0,0] = np.sum(m*(y**2 + z**2))
        self._inertia_tensor[1,1] = np.sum(m*(x**2 + z**2))
        self._inertia_tensor[2,2] = np.sum(m*(x**2 + y**2))
        self._inertia_tensor[0,1] = -np.sum(m*x*y)
        self._inertia_tensor[1,0] = self._inertia_tensor[0,1]
        self._inertia_tensor[0,2] = -np.sum(m*x*z)
        self._inertia_tensor[2,0] = self._inertia_tensor[0,2]
        self._inertia_tensor[1,2] = -np.sum(m*y*z)
        self._inertia_tensor[2,1] = self._inertia_tensor[1,2]


    def _compute_rotational_partition_function(self):
        """
        Computes the rotational temperature, the rotational constant and the rotational partition function of the system.
    
        """
        self._eigenvalues_I_SI = self._eigenvalues_I*PhysicalConstants["u"]*PhysicalConstants["A"]**2
        self._rotational_temperature = (PhysicalConstants["h"]**2/(8*np.pi**2*PhysicalConstants["kB"]))/self._eigenvalues_I_SI

        self._rotational_constant = (((PhysicalConstants["hbar"]**2)/2)/self._eigenvalues_I_SI)/PhysicalConstants["h"]*PhysicalConstants["HztoGHz"]

        self._rotational_temperature_xyz = self._rotational_temperature[0]*self._rotational_temperature[1]*self._rotational_temperature[2]
       
        self._rotational_partition_function = (np.pi**(1/2)/self._system.rotational_symmetry_number())*(self._temperature**(3/2) / (np.power(self._rotational_temperature_xyz,1/2)))
        
    def _compute_rotational_entropy(self):
        """
        Computes the rotational entropy of the system.
        """
        self._rotational_entropy = PhysicalConstants["R"]*(np.log(self._rotational_partition_function)+3/2)/PhysicalConstants["cal"]

    def _compute_rotational_energy(self):
        """
        Computes the rotational energy of the system.
        """
        self._rotational_energy = ((3/2)*PhysicalConstants["R"]*self._temperature)/PhysicalConstants["cal"]
       

    def _compute_rotational_heat_capacity(self):
        """
        Computes the rotational heat capacity of the system.
        """
        self._rotational_heat_capacity = (3/2)*PhysicalConstants["R"]/PhysicalConstants["cal"]
    
    def _rotational_contribution(self):
        self._relocate_to_cm()
        self._compute_inertia_tensor()
        self._eigenvalues_I, self._eigenvectors_I = np.linalg.eig(self._inertia_tensor) 
        self._compute_rotational_partition_function()
        self._compute_rotational_entropy()
        self._compute_rotational_energy()
        self._compute_rotational_heat_capacity()

        # print all rotational properties
        print("########################################")
        print("########################################")
        print("Rotational contribution:")
        print("Inertia tensor: ", self._inertia_tensor)
        print("Eigenvalues of inertia tensor: ", self._eigenvalues_I," in amu * Angstrom^2", " or\n", self._eigenvalues_I_SI, " in kg * m^2")
        print("Eigenvectors of inertia tensor: ", self._eigenvectors_I," in Angstrom")
        print("########################################")
        print("Rotational temperature: ", self._rotational_temperature, "in K")
        print("Rotational temperature xyz: ", self._rotational_temperature_xyz, "in K")
        print("Rotational constant: ", self._rotational_constant, "in GHz")
        print("Rotational partition function: ", self._rotational_partition_function)
        print("Rotational entropy: ", self._rotational_entropy, "in cal/(mol*K)", " or\n", self._rotational_entropy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H/T per particel")
        print("Rotational energy: ", self._rotational_energy, "in cal/mol", " or\n", self._rotational_energy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("Rotational heat capacity: ", self._rotational_heat_capacity, "in cal/(mol*K)", " or\n", self._rotational_heat_capacity*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("\n\n\n")
   

    def _compute_vibrational_partition_function(self):
        """
        Computes the vibrational temperature and the vibrational partition function of the system.
        """

        self._vib_temp_K = PhysicalConstants["h"]*self._system.real_vibrational_frequencies*PhysicalConstants["c"]*10**2/(PhysicalConstants["kB"])

        self._vibrational_temperature = np.sum(self._vib_temp_K/2)

        self._vibrational_partition_function = np.prod(np.exp(-self._vib_temp_K/(2*self._temperature))/(1-np.exp(-self._vib_temp_K/self._temperature)))

 
    def _compute_vibrational_entropy(self):
        """
        Computes the vibrational entropy of the system.
        """
    
        self._vibrational_entropy = np.subtract(np.divide((self._vib_temp_K/self._temperature),(np.exp(self._vib_temp_K/self._temperature)-1)), np.log(1-np.exp(-self._vib_temp_K/self._temperature)))

        self._vibrational_entropy = PhysicalConstants["R"]*np.sum(self._vibrational_entropy)/PhysicalConstants["cal"]

    def _compute_vibrational_energy(self):
        """
        Computes the vibrational energy of the system and the zero point energy correction
        """
        self._zpecorr = PhysicalConstants["R"]*(np.sum(self._vib_temp_K/2))/PhysicalConstants["cal"]

        self._vibrational_energy = PhysicalConstants["R"]*(np.sum(np.multiply(self._vib_temp_K,(1/2+1/(np.exp(self._vib_temp_K/self._temperature)-1)))))/PhysicalConstants["cal"]


        self._EZP = PhysicalConstants["R"]*(np.sum(self._vib_temp_K/2))/PhysicalConstants["cal"]

    def _compute_vibrational_heat_capacity(self):
        """
        Computes the vibrational heat capacity of the system.
        """
        self._vibrational_heat_capacity = PhysicalConstants["R"]*np.sum(np.exp(-self._vib_temp_K/self._temperature)*((self._vib_temp_K/self._temperature)/(np.exp(-self._vib_temp_K/self._temperature)-1))**2)/PhysicalConstants["cal"]

    def _vibrational_contribution(self):
        self._compute_vibrational_partition_function()
        self._compute_vibrational_entropy()
        self._compute_vibrational_energy()
        self._compute_vibrational_heat_capacity()

        # print all vibrational properties
        print("########################################")
        print("########################################")
        print("Vibrational temperature: ", self._vibrational_temperature, "in K")
        print("Vibrational partition function: ", self._vibrational_partition_function)
        print("Vibrational entropy: ", self._vibrational_entropy, "in cal/(mol*K)", " or\n", self._vibrational_entropy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("Vibrational energy: ", self._vibrational_energy, "in cal", " or\n", self._vibrational_energy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle", " or\n", self._vibrational_energy/PhysicalConstants["cal"]/1000, " in kcal")
        print("Vibrational heat capacity: ", self._vibrational_heat_capacity, "in cal/(mol*K)", " or\n", self._vibrational_heat_capacity*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H/T per particle")
        print("######################################### ")
        print("Zero point energy correction: ", self._zpecorr, "in cal", " or\n", self._zpecorr*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle", " or\n", self._zpecorr/PhysicalConstants["cal"]/1000, " in kcal")
        print("Zero point vibrational energy: ", self._EZP, "in cal/mol", " or\n", self._EZP*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle", " in cal/mol", " or\n", self._EZP/PhysicalConstants["cal"]/1000, " in kcal/mol")
        print("\n\n\n")

    def _compute_electronic_partition_function(self):
        """
        Computes the electronic partition function of the system.
        """
    
        self._electronic_partition_function = 2 * self._system._spin + 1

    def _compute_electronic_entropy(self):
        """
        Computes the electronic entropy of the system.
        """
        self._electronic_entropy = PhysicalConstants["R"]*np.log(self._electronic_partition_function)/PhysicalConstants["cal"]
    
    def _compute_electronic_energy(self):
        """
        Computes the electronic energy of the system.
        """
        self._electronic_energy = 0/PhysicalConstants["cal"]
    
    def _compute_electronic_heat_capacity(self):
        """
        Computes the electronic heat capacity of the system.
        """
        self._electronic_heat_capacity = 0/PhysicalConstants["cal"]
    
    def _electronic_contribution(self):
        self._compute_electronic_partition_function()
        self._compute_electronic_entropy()
        self._compute_electronic_energy()
        self._compute_electronic_heat_capacity()

        # print all electronic properties
        print("########################################")
        print("########################################")
        print("Electronic partition function: ", self._electronic_partition_function)
        print("Electronic entropy: ", self._electronic_entropy, "in cal/(mol*K)", " or\n", self._electronic_entropy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H/T per particle")
        print("Electronic energy: ", self._electronic_energy, "in cal/mol", " or\n", self._electronic_energy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("Electronic heat capacity: ", self._electronic_heat_capacity, "in cal/(mol*K)", " or\n", self._electronic_heat_capacity*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("\n\n\n")

    def _compute_translational_partition_function(self):
        molecular_mass = np.sum(self._system.atomic_masses())*PhysicalConstants["u"]

        self._translational_partition_function = np.power((2*np.pi*molecular_mass*PhysicalConstants["kB"]*self._temperature)/(PhysicalConstants["h"]**2),(3/2))*(PhysicalConstants["kB"]*self._temperature/self._pressure)
    
    def _compute_translational_entropy(self):
        self._translational_entropy = PhysicalConstants["R"]*(np.log(self._translational_partition_function)+5/2)/PhysicalConstants["cal"]
    
    def _compute_translational_energy(self):
        self._translational_energy = 3/2*PhysicalConstants["R"]*self._temperature/PhysicalConstants["cal"]

    def _compute_translational_heat_capacity(self):
        self._translational_heatcapacity = 3/2*PhysicalConstants["R"]/PhysicalConstants["cal"]
    
    def _translational_contribution(self):
        self._compute_translational_partition_function()
        self._compute_translational_entropy()
        self._compute_translational_energy()
        self._compute_translational_heat_capacity()

        # print all translational properties
        print("########################################")
        print("########################################")
        print("Translational partition function: ", self._translational_partition_function)
        print("Translational entropy: ", self._translational_entropy, "in cal/(mol*K)", " or\n", self._translational_entropy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("Translational energy: ", self._translational_energy, "in cal", " or\n", self._translational_energy*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle", " or\n", self._translational_energy/PhysicalConstants["cal"]/1000, " in kcal")
        print("Translational heat capacity: ", self._translational_heatcapacity, "in cal/(mol*K)", " or\n", self._translational_heatcapacity*PhysicalConstants["cal"]/PhysicalConstants["H"]/PhysicalConstants["Na"], " in H per particle")
        print("\n\n\n")


    def _compute_thermochemical_properties(self):
        """
        Computes the thermochemical properties of the system.
        """
        self._total_energy  = self._rotational_energy + self._vibrational_energy + self._electronic_energy + self._translational_energy
        self._total_energy_kcal = self._total_energy/1000
        self._total_energy_Hartree_per_mol = self._total_energy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]

       
        self._total_entropy = self._rotational_entropy + self._vibrational_entropy + self._electronic_entropy + self._translational_entropy
        self._total_entropy_Hartree_per_mol = self._total_entropy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]


        self._total_enthalpy = self._total_energy + (PhysicalConstants["R"]*self._temperature)/PhysicalConstants["cal"]
        self._total_enthalpy_kcal = self._total_enthalpy/1000
        self._total_enthalpy_Hartree_per_mol = self._total_enthalpy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]


        self._total_gibbs_free_energy = self._total_energy - self._total_entropy * self._temperature
        self._total_gibbs_free_energy_kcal = self._total_gibbs_free_energy/1000
        self._total_gibbs_free_energy_Hartree_per_mol = self._total_gibbs_free_energy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]
        
        
        self._total_heatcapacity = self._rotational_heat_capacity + self._vibrational_heat_capacity + self._electronic_heat_capacity + self._translational_heatcapacity

    def _summary(self):
        self._EeZPE = self._system.electronic_energy+ self._zpecorr/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]
        self._EeEtot = self._system.electronic_energy + self._total_energy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]
        self._EeHtot = self._system.electronic_energy + self._total_enthalpy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]
        self._EeGtot = self._system.electronic_energy + self._total_gibbs_free_energy/PhysicalConstants["JinHartree"]/PhysicalConstants["Na"]*PhysicalConstants["cal"]

        
        

        print("########################################")
        print("########################################")
        print("#################Summary:###############")
        print("Total energy: ", self._total_energy, " in cal", " or\n", self._total_energy_kcal, " in kcal", " or\n", self._total_energy_Hartree_per_mol, " in Hartree per mol")
        print("Total entropy: ", self._total_entropy, "in cal/(mol*K)", " or\n", self._total_entropy_Hartree_per_mol, " in Hartree per mol")
        print("Total enthalpy: ", self._total_enthalpy, "in cal", " or\n", self._total_enthalpy_kcal, " in kcal", " or\n", self._total_enthalpy_Hartree_per_mol, " in Hartree per mol")
        print("Total Gibbs free energy: ", self._total_gibbs_free_energy, " in cal", " or\n", self._total_gibbs_free_energy_kcal, " in kcal", " or\n", self._total_gibbs_free_energy_Hartree_per_mol, " in Hartree per mol")
        print("Total heat capacity: ", self._total_heatcapacity, " in cal/(mol*K)")
        print("########################################")
        print("########################################")
        print("########################################")
        print("Sum of electronic energy and zero point energy correction: \n", self._EeZPE, "Hartree per particle")
        print("Sum of electronic energy and total energy: \n", self._EeEtot, "Hartree per particle")
        print("Sum of electronic energy and total enthalpy: \n", self._EeHtot, "Hartree per particle")
        print("Sum of electronic energy and total Gibbs free energy: \n", self._EeGtot, "Hartree per particle")
        print("########################################")
        print("########################################")
        print("########################################")
