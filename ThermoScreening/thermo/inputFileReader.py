class InputFileReader:
    """
    Class for Reading the Input File. It reads the input file and checks if all required keys are set.

    Parameters
    ----------
    input_file : str
        The input file.

    """

    # required keys
    required_keys = [
        "coord_file",
        "temperature",
        "pressure",
        "engine",
        "vibrational_file",
        "energy",
    ]

    def __init__(
        self,
        input_file: str | None = None,
    ) -> None:
        """
        Initializes the InputFileReader class.

        Parameters
        ----------
        _input_file : str, optional, default=None
            The input file.

        Raises
        ------
        ValueError
            If the input file is not given.

        Returns
        -------
        None
        """
        if input_file is None:
            raise ValueError(
                "The input file has to be given to initialize the InputFileReader."
            )
        self._input_file = input_file
        self._read()
        self._check()

        return None

    def _read(self):
        """
        Reads the input file and parses it.
        It also sets the raw_input_file and the dictionary.

        Returns
        -------
        None
        """
        self._raw_input_file = open(self._input_file, "r").readlines()
        self._dictionary = {}
        for line in self._raw_input_file:
            if line[0] != "#":
                key, value = line.split(" = ")
                self._dictionary[key.strip()] = value.strip()
                print(key.strip(), " = ", value.strip())

        return None

    def _check(self):
        """
        Checks if all required keys are set and if all keys are known.

        Raises
        ------
        ValueError
            If a required key is not set or if an unknown key is set.
        
        Returns
        -------
        None
        """
        self._check_required_keys()
        self._check_known_keys()
        
        return None

    def _check_required_keys(self):
        """
        Checks if all required keys are set.

        Raises
        ------
        ValueError
            If a required key is not set.

        Returns
        -------
        None
        """
        for key in self.required_keys:
            if key not in self._dictionary.keys():
                raise ValueError("The key {} is not set in the input file.".format(key))
            
        return None

    def _check_known_keys(self):
        """
        Checks if all keys are known. 

        Raises
        ------
        ValueError
            If an unknown key is set.

        Returns
        -------
        None
        """
        for key in self._dictionary.keys():
            if key not in self.required_keys:
                raise ValueError("The key {} is not known.".format(key))

        return None