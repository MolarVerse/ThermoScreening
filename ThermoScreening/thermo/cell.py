import numpy as np

"""
Dummy class for cell
"""
def cell(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
        """
        calculate cell vectors from a, b, c, alpha, beta, gamma#

        Parameters
        ----------
        a : float
            The length of the first cell vector.
        b : float
            The length of the second cell vector.
        c : float
            The length of the third cell vector.s
        alpha : float
            The angle between the second and third cell vectors.
        beta : float
            The angle between the first and third cell vectors.
        gamma : float
            The angle between the first and second cell vectors.

        """
        volume = a*b*c*np.sqrt(1+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)-np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2)
        cell_vector = np.zeros((3,3))
        cell_vector[0,0] = a
        cell_vector[1,0] = 0.0
        cell_vector[2,0] = 0.0
        cell_vector[0,1] = b*np.cos(gamma)
        cell_vector[1,1] = b*np.sin(gamma)
        cell_vector[2,1] = 0.0
        cell_vector[0,2] = c*np.cos(beta)
        cell_vector[1,2] = c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
        cell_vector[2,2] = volume/(a*b*np.sin(gamma))
        return cell_vector
def cell_parameters_calc(cell_vector: np.ndarray) -> np.ndarray:  
    """
    calculate cell parameters from cell vectors

    Parameters
    ----------
    cell_vector : np.ndarray
        The cell vectors.

    return
    ------
    np.ndarray
        The cell parameters.

    """
    print(cell_vector)
    a = np.linalg.norm(cell_vector[:,0])
    b = np.linalg.norm(cell_vector[:,1])
    c = np.linalg.norm(cell_vector[:,2])
    alpha = np.arccos(np.dot(cell_vector[:,1],cell_vector[:,2])/(b*c))
    beta = np.arccos(np.dot(cell_vector[:,0],cell_vector[:,2])/(a*c))
    gamma = np.arccos(np.dot(cell_vector[:,0],cell_vector[:,1])/(a*b))
    return np.array([a, b, c, alpha, beta, gamma])
     
     

class Cell:
    """
    Dummy class for cell
    """
    def __init__(self,
                cell_vectors : np.ndarray | None = None,
                a : float | None = None,
                b : float | None = None,
                c : float | None = None,
                alpha : float | None = None,
                beta : float | None = None,
                gamma : float | None = None,
                volume : float | None = None
                ) -> None:
        """
        Initializes the Cell class.

        Parameters
        ----------
        cell_vectors : np.ndarray, optional, default=None
           The cell vectors.
        a : float, optional, default=None
           The length of the first cell vector.
        b : float, optional, default=None
           The length of the second cell vector.
        c : float, optional, default=None
           The length of the third cell vector.
        alpha : float, optional, default=None
           The angle between the second and third cell vectors.
        beta : float, optional, default=None
           The angle between the first and third cell vectors.
        gamma : float, optional, default=None
           The angle between the first and second cell vectors.
        volume : float, optional, default=None
            The volume of the cell.
        """
   
        self._a = a
        self._b = b
        self._c = c
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma
        self._volume = volume
        self._cell_parameters = np.array([a, b, c, alpha, beta, gamma])
        self._cell_vectors = cell_vectors
        self._cell_inverse = None
       
        
        if self._cell_parameters is None:
            self._cell_parameters = cell_parameters_calc(cell_vectors)
        else:
            self._cell_vectors = cell(self._cell_parameters[0], self._cell_parameters[1], self._cell_parameters[2], self._cell_parameters[3], self._cell_parameters[4], self._cell_parameters[5])
        
        if self._volume is None:
            self._volume = self._cell_parameters[0]*self._cell_parameters[1]*self._cell_parameters[2]*np.sqrt(1+2*np.cos(self._cell_parameters[3])*np.cos(self._cell_parameters[4])*np.cos(self._cell_parameters[5])-np.cos(self._cell_parameters[3])**2-np.cos(self._cell_parameters[4])**2-np.cos(self._cell_parameters[5])**2)
        
        

        if self._cell_inverse is None:
             self._cell_inverse = np.linalg.inv(self._cell_vectors)

    @property
    def cell_vectors(self) -> np.ndarray:
        """
        Cell vectors.

        Returns
        -------
        np.ndarray
            The cell vectors.
        """
        return self._cell_vectors
    
    @property
    def cell_inverse(self) -> np.ndarray:
        """
        Inverse of the cell vectors.

        Returns
        -------
        np.ndarray
            The inverse of the cell vectors.
        """
        return self._cell_inverse
    
    @property
    def cell_parameters(self) -> np.ndarray:
        """
        Cell parameters.

        Returns
        -------
        np.ndarray
            The cell parameters.
        """
        return self._cell_parameters
    
    @property
    def a(self) -> float:
        """
        Length of the first cell vector.

        Returns
        -------
        float
            The length of the first cell vector.
        """
        return self._a
    
    @property
    def b(self) -> float:
        """
        Length of the second cell vector.

        Returns
        -------
        float
            The length of the second cell vector.
        """
        return self._b
    
    @property
    def c(self) -> float:
        """
        Length of the third cell vector.

        Returns
        -------
        float
            The length of the third cell vector.
        """
        return self._c
    
    @property
    def alpha(self) -> float:
        """
        Angle between the second and third cell vectors.

        Returns
        -------
        float
            The angle between the second and third cell vectors.
        """
        return self._alpha
    
    @property
    def beta(self) -> float:
        """
        Angle between the first and third cell vectors.

        Returns
        -------
        float
            The angle between the first and third cell vectors.
        """
        return self._beta
    
    @property
    def gamma(self) -> float:
        """
        Angle between the first and second cell vectors.

        Returns
        -------
        float
            The angle between the first and second cell vectors.
        """
        return self._gamma
    

    @property
    def volume(self) -> float:
        """
        Volume of the cell.

        Returns
        -------
        float
            The volume of the cell.
        """
        return self._volume
    
    def __repr__(self) -> str:
        """
        Returns a string representation of the Cell class.

        Returns
        -------
        str
            The string representation of the Cell class.
        """
        return f"Cell(a={self._a}, b={self._b}, c={self._c}, alpha={self._alpha}, beta={self._beta}, gamma={self._gamma})"
    
