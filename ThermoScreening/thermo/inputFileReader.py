import logging

from ThermoScreening.exceptions import TSValueError
from ThermoScreening.utils.custom_logging import setup_logger
from ThermoScreening import __package_name__

class InputFileReader:
    """
    Class for Reading the Input File. It reads the input file and checks if all required keys are set.

    Parameters
    ----------
    input_file : str
        The input file.

    """
    
    logger = logging.getLogger(__package_name__).getChild(__qualname__)
    logger = setup_logger(logger)

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
        TSValueError
            If the input file is not given.

        Returns
        -------
        None
        """
        if input_file is None:
            self.logger.error(
                "The input file has to be given to initialize the InputFileReader.",
                exception=TSValueError,
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
        TSValueError
            If a required key is not set.

        Returns
        -------
        None
        """
        for key in self.required_keys:
            if key not in self._dictionary.keys():
                self.logger.error(
                    "The key {} is not set in the input file.".format(key),
                    exception=TSValueError
                )
            
        return None

    def _check_known_keys(self):
        """
        Checks if all keys are known. 

        Raises
        ------
        TSValueError
            If an unknown key is set.

        Returns
        -------
        None
        """
        for key in self._dictionary.keys():
            if key not in self.required_keys:
                self.logger.error(
                    "The key {} is not known.".format(key),
                    exception=TSValueError
                )

        return None