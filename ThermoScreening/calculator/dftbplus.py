from ase.calculators.dftb import Dftb


class Geoopt(Dftb):
    """
    Custom DFTB+ calculator. It is a subclass of ase.calculators.dftb.Dftb.

    Parameters:
    -----------
    atoms : ase.Atoms
        Atoms object.
    label : str
        Label for the calculation. Default is 'geo_opt'.
    charge : int
        Charge of the system. Default is 0.
    slako_dir : str
        Path to the Slater-Koster files. If None, it will look for the DFTB_PREFIX environment variable.

    Other Parameters:
    -----------------
    **kwargs : dict
        Additional keyword arguments to pass to the Dftb class.

    Methods:
    --------
    run():
        Run the geometry optimization.
    """

    def __init__(self, atoms, label='geo_opt', charge=0, slako_dir=None, **kwargs):
        super().__init__(
            atoms=atoms,
            label=label,
            slako_dir=slako_dir,
            Hamiltonian_Charge=charge,
            Driver_='LBFGS', 
            Driver_MaxForceComponent='1.0e-5',
            Driver_OutputPrefix='geo_opt',
            **kwargs
        )
        self.atoms = atoms
        self.label = label
        self.charge = charge
        self.slako_dir = slako_dir
        self.kwargs = kwargs

        return None

    def run(self):
        """
        Run the geometry optimization.

        Returns:
        --------
        ase.Atoms
            Atoms object with the final geometry.
        """
        self.atoms.set_calculator(self)
        self.atoms.get_potential_energy()

        return self.atoms
