import numpy as np
from scipy.io import FortranFile

def main():
    endial='<'
    xgrid = 11
    ygrid = 11
    filename = 'out.dat'
    #filename = 'fort.125'
    f = FortranFile(filename, 'r')
    t = f.read_record('f8')[0]
    U = f.read_record('({},{})f8'.format(ygrid, xgrid)).T
    V = f.read_record('({},{})f8'.format(ygrid, xgrid)).T
    X = f.read_record('f8')
    Y = f.read_record('f8')
    print(t)
    print(U)
    print(V)
    print(X)
    print(Y)
    f.close()

if __name__ == '__main__':
    main()
