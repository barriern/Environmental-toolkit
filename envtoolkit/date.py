
import numpy as np

def leapy(year):

    """ Returns whether a year is leap or not 

    :param int year: Year
    :return: True if leap, else False
    :rtype: bool
    
    """
    
    if (year % 4 == 0) & (year % 100 == 0):
        return True

    else:
        if year % 400 == 0:
            return True
        else:
            return False




if __name__ == "__main__":

    print leapy(2000)
    x = np.arange(1990, 2001)
    x = np.tile(x, (4, 1))
    print x.shape
    print leapy(x)

