import numpy as np
from ..utils.physicalConstants import PhysicalConstants
from .system import System


class Thermo:
    """
    A class for calculating the thermochemistral properties e.g. Gibbs free energy,
    enthalpy, entropy, heat capacity, etc. of a system. The class needs as input the
    temperature, pressure, system information (coordinates, degree of freedom,
    electronic spin, symmetry number, electronic energy and vibrational frequencies).

    Attributes
    ----------
    temperature : float
        The temperature of the system in Kelvin.
    pressure : float
        The pressure of the system in bar.
    system : System
        The system information of the system.
    engine : str
        The engine used for the calculation to compute the thermochemical properties
        with correct units.

    Methods
    -------
    run()
        Calculates the thermochemical properties of the system.
    """

    def __init__(
        self, temperature: float, pressure: float, system: System, engine: str
    ):
        """
        Initializes the Thermo class with the temperature, pressure, system information

        Parameters
        ----------
        temperature : float
            The temperature of the system in Kelvin.
        pressure : float
            The pressure of the system in bar.
        system : System
            The system information of the system.
        engine : str
            The engine used for the calculation to compute the thermochemical
            properties with correct units.

        Raises
        ------
        ValueError
            If the pressure is not given.
            If the temperature is not given.
            If the system is not given.
            If the engine is not supported.
            If the temperature is negative.
            If the pressure is negative.

        Returns
        -------
        None
        """
        if pressure == None:
            raise ValueError("The pressure is not given.")
        if temperature == None:
            raise ValueError("The temperature is not given.")
        if system == None:
            raise ValueError("The system is not given.")
        if engine == None:
            raise ValueError("The engine is not given.")

        self._temperature = temperature
        self._pressure = pressure
        self._system = system
        self._engine = engine

        if self._engine != "dftb+":
            raise ValueError("The engine is not supported.")
        if self._temperature < 0:
            raise ValueError("The temperature is negative.")
        if self._pressure < 0:
            raise ValueError("The pressure is negative.")

        self._coord = self._system.coord()
        self._atomic_masses = self._system.atomic_masses()

        return None

    def run(self):
        """
        Calculates the thermochemical properties of the system.

        Returns
        -------
        None
        """
        # self._transform_units()
        self._rotational_contribution()
        self._vibrational_contribution()
        self._electronic_contribution()
        self._translational_contribution()

        self._compute_thermochemical_properties()
        self._summary()

        return None

    def _transform_units(self):
        """
        Transforms the units of the system to the correct units.

        Returns
        -------
        None
        """
        if self._engine == "DFTB+":
            # the originial units of DFTB+ are Hartree, Angstrom, amu
            # the units are transformed to J, m, kg
            self._coord = self._system.coord()[:] * PhysicalConstants["A"]
            self._atomic_masses = (
                self._system.atomic_masses()[:] * PhysicalConstants["u"]
            )
            self._system.real_vibrational_frequencies[:] = (
                self._system.real_vibrational_frequencies[:]
                * PhysicalConstants["HztoGHz"]
            )
            self._system.electronic_energy = (
                self._system.electronic_energy
                * PhysicalConstants["H"]
                * PhysicalConstants["N_A"]
            )
            self._system.real_vibrational_frequencies[:] = (
                self._system.real_vibrational_frequencies[:]
                * PhysicalConstants["HztoGHz"]
            )

        else:
            raise ValueError("The engine is not supported.")

        return None

    def _relocate_to_cm(self):
        """
        Relocates the system to the center of mass.

        Returns
        -------
        None
        """

        # ? Resolving the issue with the center of mass calculation
        coord = self._system.coord()
        coord -= self._system.center_of_mass

        return None

    def _compute_inertia_tensor(self):
        """
        Computes the inertia tensor of the system.

        Returns
        -------
        None
        """
        coord = self._system.coord()
        x = coord[:, 0]
        y = coord[:, 1]
        z = coord[:, 2]

        m = self._system.atomic_masses()

        self._inertia_tensor = np.zeros((3, 3))

        self._inertia_tensor[0, 0] = np.sum(m * (y**2 + z**2))
        self._inertia_tensor[1, 1] = np.sum(m * (x**2 + z**2))
        self._inertia_tensor[2, 2] = np.sum(m * (x**2 + y**2))
        self._inertia_tensor[0, 1] = -np.sum(m * x * y)
        self._inertia_tensor[1, 0] = self._inertia_tensor[0, 1]
        self._inertia_tensor[0, 2] = -np.sum(m * x * z)
        self._inertia_tensor[2, 0] = self._inertia_tensor[0, 2]
        self._inertia_tensor[1, 2] = -np.sum(m * y * z)
        self._inertia_tensor[2, 1] = self._inertia_tensor[1, 2]

        return None

    def _compute_rotational_partition_function(self):
        """
        Computes the rotational temperature, the rotational constant and the rotational partition function of the system.

        Returns
        -------
        None
        """

        # TODO: resolve with naming pylint issue 
        self._eigenvalues_I_SI = (
            self._eigenvalues_I * PhysicalConstants["u"] * PhysicalConstants["A"] ** 2
        )
        self._rotational_temperature = (
            PhysicalConstants["h"] ** 2 / (8 * np.pi**2 * PhysicalConstants["kB"])
        ) / self._eigenvalues_I_SI

        self._rotational_constant = (
            (((PhysicalConstants["hbar"] ** 2) / 2) / self._eigenvalues_I_SI)
            / PhysicalConstants["h"]
            * PhysicalConstants["HztoGHz"]
        )

        self._rotational_temperature_xyz = (
            self._rotational_temperature[0]
            * self._rotational_temperature[1]
            * self._rotational_temperature[2]
        )

        self._rotational_partition_function = (
            np.pi ** (1 / 2) / self._system.rotational_symmetry_number()
        ) * (
            self._temperature ** (3 / 2)
            / (np.power(self._rotational_temperature_xyz, 1 / 2))
        )

        return None

    def _compute_rotational_entropy(self):
        """
        Computes the rotational entropy of the system.

        Returns
        -------
        None
        """
        self._rotational_entropy = (
            PhysicalConstants["R"]
            * (np.log(self._rotational_partition_function) + 3 / 2)
            / PhysicalConstants["cal"]
        )

        return None

    def _compute_rotational_energy(self):
        """
        Computes the rotational energy of the system.

        Returns
        -------
        None
        """
        self._rotational_energy = (
            (3 / 2) * PhysicalConstants["R"] * self._temperature
        ) / PhysicalConstants["cal"]

    def _compute_rotational_heat_capacity(self):
        """
        Computes the rotational heat capacity of the system.
        """
        self._rotational_heat_capacity = (
            (3 / 2) * PhysicalConstants["R"] / PhysicalConstants["cal"]
        )

        return None

    def _rotational_contribution(self):
        """
        Computes the rotational contribution of the system.

        Returns
        -------
        None
        """

        self._relocate_to_cm()
        self._compute_inertia_tensor()
        self._eigenvalues_I, self._eigenvectors_I = np.linalg.eig(self._inertia_tensor)
        self._compute_rotational_partition_function()
        self._compute_rotational_entropy()
        self._compute_rotational_energy()
        self._compute_rotational_heat_capacity()

        return None

    def _compute_vibrational_partition_function(self):
        """
        Computes the vibrational temperature and the vibrational partition function of the system.

        Returns
        -------
        None
        """
        self._vib_temp_K = (
            PhysicalConstants["h"]
            * self._system.real_vibrational_frequencies
            * PhysicalConstants["c"]
            * 10**2
            / (PhysicalConstants["kB"])
        )

        self._vibrational_temperature = np.sum(self._vib_temp_K / 2)

        self._vibrational_partition_function = np.prod(
            np.exp(-self._vib_temp_K / (2 * self._temperature))
            / (1 - np.exp(-self._vib_temp_K / self._temperature))
        )

        return None

    def _compute_vibrational_entropy(self):
        """
        Computes the vibrational entropy of the system.

        Returns
        -------
        None
        """

        self._vibrational_entropy = np.subtract(
            np.divide(
                (self._vib_temp_K / self._temperature),
                (np.exp(self._vib_temp_K / self._temperature) - 1),
            ),
            np.log(1 - np.exp(-self._vib_temp_K / self._temperature)),
        )

        self._vibrational_entropy = (
            PhysicalConstants["R"]
            * np.sum(self._vibrational_entropy)
            / PhysicalConstants["cal"]
        )

        return None

    def _compute_vibrational_energy(self):
        """
        Computes the vibrational energy of the system and the zero point energy correction.

        Returns
        -------
        None
        """

        self._zpecorr = (
            PhysicalConstants["R"]
            * (np.sum(self._vib_temp_K / 2))
            / PhysicalConstants["cal"]
        )

        self._vibrational_energy = (
            PhysicalConstants["R"]
            * (
                np.sum(
                    np.multiply(
                        self._vib_temp_K,
                        (
                            1 / 2
                            + 1 / (np.exp(self._vib_temp_K / self._temperature) - 1)
                        ),
                    )
                )
            )
            / PhysicalConstants["cal"]
        )

        self._EZP = (
            PhysicalConstants["R"]
            * (np.sum(self._vib_temp_K / 2))
            / PhysicalConstants["cal"]
        )

        return None

    def _compute_vibrational_heat_capacity(self):
        """
        Computes the vibrational heat capacity of the system.

        Returns
        -------
        None
        """
        self._vibrational_heat_capacity = (
            PhysicalConstants["R"]
            * np.sum(
                np.exp(-self._vib_temp_K / self._temperature)
                * (
                    (self._vib_temp_K / self._temperature)
                    / (np.exp(-self._vib_temp_K / self._temperature) - 1)
                )
                ** 2
            )
            / PhysicalConstants["cal"]
        )
        return None

    def _vibrational_contribution(self):
        """
        Computes the vibrational contribution of the system.

        Returns
        -------
        None
        """

        self._compute_vibrational_partition_function()
        self._compute_vibrational_entropy()
        self._compute_vibrational_energy()
        self._compute_vibrational_heat_capacity()

    def _compute_electronic_partition_function(self):
        """
        Computes the electronic partition function of the system.

        Returns
        -------
        None
        """

        self._electronic_partition_function = 2 * self._system._spin + 1

        return None

    def _compute_electronic_entropy(self):
        """
        Computes the electronic entropy of the system.

        Returns
        -------
        None
        """

        self._electronic_entropy = (
            PhysicalConstants["R"]
            * np.log(self._electronic_partition_function)
            / PhysicalConstants["cal"]
        )
        return None

    def _compute_electronic_energy(self):
        """
        Computes the electronic energy of the system.

        Returns
        -------
        None
        """

        self._electronic_energy = 0 / PhysicalConstants["cal"]

        return None

    def _compute_electronic_heat_capacity(self):
        """
        Computes the electronic heat capacity of the system.

        Returns
        -------
        None
        """

        self._electronic_heat_capacity = 0 / PhysicalConstants["cal"]

        return None

    def _electronic_contribution(self):
        """
        Computes the electronic contribution of the system.

        Returns
        -------
        None
        """

        self._compute_electronic_partition_function()
        self._compute_electronic_entropy()
        self._compute_electronic_energy()
        self._compute_electronic_heat_capacity()

        return None

    def _compute_translational_partition_function(self):
        """
        Computes the translational partition function of the system.

        Returns
        -------
        None
        """

        molecular_mass = np.sum(self._system.atomic_masses()) * PhysicalConstants["u"]

        self._translational_partition_function = np.power(
            (2 * np.pi * molecular_mass * PhysicalConstants["kB"] * self._temperature)
            / (PhysicalConstants["h"] ** 2),
            (3 / 2),
        ) * (PhysicalConstants["kB"] * self._temperature / self._pressure)

        return None

    def _compute_translational_entropy(self):
        """
        Computes the translational entropy of the system.

        Returns
        -------
        None
        """

        self._translational_entropy = (
            PhysicalConstants["R"]
            * (np.log(self._translational_partition_function) + 5 / 2)
            / PhysicalConstants["cal"]
        )

        return None

    def _compute_translational_energy(self):
        """
        Computes the translational energy of the system.

        Returns
        -------
        None
        """

        self._translational_energy = (
            3
            / 2
            * PhysicalConstants["R"]
            * self._temperature
            / PhysicalConstants["cal"]
        )

        return None

    def _compute_translational_heat_capacity(self):
        """
        Computes the translational heat capacity of the system.

        Returns
        -------
        None
        """

        self._translational_heatcapacity = (
            3 / 2 * PhysicalConstants["R"] / PhysicalConstants["cal"]
        )

        return None

    def _translational_contribution(self):
        """
        Computes the translational contribution of the system.

        Returns
        -------
        None
        """
        self._compute_translational_partition_function()
        self._compute_translational_entropy()
        self._compute_translational_energy()
        self._compute_translational_heat_capacity()

        return None

    def _compute_thermochemical_properties(self):
        """
        Computes the thermochemical properties of the system.

        Returns
        -------
        None
        """
        self._total_energy = (
            self._rotational_energy
            + self._vibrational_energy
            + self._electronic_energy
            + self._translational_energy
        )
        self._total_energy_kcal = self._total_energy / 1000
        self._total_energy_Hartree_per_mol = (
            self._total_energy
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )

        self._total_entropy = (
            self._rotational_entropy
            + self._vibrational_entropy
            + self._electronic_entropy
            + self._translational_entropy
        )
        self._total_entropy_Hartree_per_mol = (
            self._total_entropy
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )

        self._total_enthalpy = (
            self._total_energy
            + (PhysicalConstants["R"] * self._temperature) / PhysicalConstants["cal"]
        )
        self._total_enthalpy_kcal = self._total_enthalpy / 1000
        self._total_enthalpy_Hartree_per_mol = (
            self._total_enthalpy
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )

        self._total_gibbs_free_energy = (
            self._total_energy - self._total_entropy * self._temperature
        )
        self._total_gibbs_free_energy_kcal = self._total_gibbs_free_energy / 1000
        self._total_gibbs_free_energy_Hartree_per_mol = (
            self._total_gibbs_free_energy
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )

        self._total_heatcapacity = (
            self._rotational_heat_capacity
            + self._vibrational_heat_capacity
            + self._electronic_heat_capacity
            + self._translational_heatcapacity
        )

        return None

    def _summary(self):
        """
        Summary of the thermochemical properties of the system.

        Returns
        -------
        None
        """

        self._EeZPE = (
            self._system.electronic_energy
            + self._zpecorr
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )
        self._EeEtot = (
            self._system.electronic_energy
            + self._total_energy
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )
        self._EeHtot = (
            self._system.electronic_energy
            + self._total_enthalpy
            / PhysicalConstants["H"]
            / PhysicalConstants["N_A"]
            * PhysicalConstants["cal"]
        )
        self._EeGtot = (
            self._system.electronic_energy
            + self._total_gibbs_free_energy_Hartree_per_mol
        )

        return None

    def total_energy(self, unit: str) -> float:
        """
        Returns the total energy of the system based on the unit.

        Parameters
        ----------
        unit : str
            The unit of the total energy.

        Raises
        ------
        ValueError
            If the unit is not supported.

        Returns
        -------
        float
            The total energy of the system.
        """

        if unit == "H":
            return self._total_energy_Hartree_per_mol
        elif unit == "kcal":
            return self._total_energy_kcal
        elif unit == "cal":
            return self._total_energy
        else:
            raise ValueError("The unit is not supported.")

    def total_enthalpy(self, unit: str) -> float:
        """
        Returns the total enthalpy of the system based on the unit.

        Parameters
        ----------
        unit : str
            The unit of the total enthalpy.

        Raises
        ------
        ValueError
            If the unit is not supported.

        Returns
        -------
        float
            The total enthalpy of the system.
        """

        if unit == "H":
            return self._total_enthalpy_Hartree_per_mol
        elif unit == "kcal":
            return self._total_enthalpy_kcal
        elif unit == "cal":
            return self._total_enthalpy
        else:
            raise ValueError("The unit is not supported.")

    def total_entropy(self, unit: str) -> float:
        """
        Returns the total entropy of the system based on the unit.

        Parameters
        ----------
        unit : str
            The unit of the total entropy.

        Raises
        ------
        ValueError
            If the unit is not supported.

        Returns
        -------
        float
            The total entropy of the system.
        """

        if unit == "H/T":
            return self._total_entropy_Hartree_per_mol
        elif unit == "cal/(mol*K)":
            return self._total_entropy
        else:
            raise ValueError("The unit is not supported.")

    def total_gibbs_free_energy(self, unit: str) -> float:
        """
        Returns the total Gibbs free energy of the system based on the unit.

        Parameters
        ----------
        unit : str
            The unit of the total Gibbs free energy.

        Raises
        ------
        ValueError
            If the unit is not supported.

        Returns
        -------
        float
            The total Gibbs free energy of the system.
        """

        if unit == "H":
            return self._total_gibbs_free_energy_Hartree_per_mol
        elif unit == "kcal":
            return self._total_gibbs_free_energy_kcal
        elif unit == "cal":
            return self._total_gibbs_free_energy
        else:
            raise ValueError("The unit is not supported.")

    def total_heat_capacity(self, unit: str) -> float:
        """
        Returns the total heat capacity of the system based on the unit.

        Parameters
        ----------
        unit : str
            The unit of the total heat capacity.

        Raises
        ------
        ValueError
            If the unit is not supported.

        Returns
        -------
        float
            The total heat capacity of the system.
        """

        if unit == "cal/(mol*K)":
            return self._total_heatcapacity
        elif unit == "H/T":
            return (
                self._total_heatcapacity
                * PhysicalConstants["cal"]
                / PhysicalConstants["H"]
                / PhysicalConstants["N_A"]
            )
        else:
            raise ValueError("The unit is not supported.")

    def total_EeZPE(self) -> float:
        """
        Returns the sum of the electronic energy and the zero point energy correction.

        Returns
        -------
        float
            The sum of the electronic energy and the zero point energy correction in Hartree per particle.
        """

        return self._EeZPE

    def total_EeEtot(self) -> float:
        """
        Returns the sum of the electronic energy and the total energy.

        Returns
        -------
        float
            The sum of the electronic energy and the total energy in Hartree per particle.
        """

        return self._EeEtot

    def total_EeHtot(self) -> float:
        """
        Returns the sum of the electronic energy and the total enthalpy.

        Returns
        -------
        float
            The sum of the electronic energy and the total enthalpy in Hartree per particle.
        """

        return self._EeHtot

    def total_EeGtot(self) -> float:
        """
        Returns the sum of the electronic energy and the total Gibbs free energy.

        Returns
        -------
        float
            The sum of the electronic energy and the total Gibbs free energy in Hartree per particle.
        """

        return self._EeGtot

    def electronic_energy(self) -> float:
        """
        Returns the electronic energy of the system.

        Returns
        -------
        float
            The electronic energy of the system in Hartree per particle.
        """

        return self._system.electronic_energy
