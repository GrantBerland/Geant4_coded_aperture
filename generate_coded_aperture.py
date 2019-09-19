import numpy as np

boxSize = 0.22;
dimX = 4.;
dimY = 4.;

#centeringShift = 0.2;
centeringShift = 0.22

# p, must be prime and satisfy L = 4*m + 1, for m in Z
gridSizeX = 17 
gridSizeY = gridSizeX; # (square)

def jacobi(a, n):
    assert(n > a > 0 and n%2 == 1)
    t = 1
    while a != 0:
        while a % 2 == 0:
            a /= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    else:
        return 0

def MURA_code(i, j, p):
    '''
    This function generates a Modified Uniform Redundent Aperature
    '''
    if i == 0:
        return 0;
    elif j == 0 and i != 0:
        return 1;
    elif jacobi(i, p)*jacobi(j, p) == 1: 
        return 1;
    else:
        return 0;

def MURA_decoding_matrix(i,j, block):
    if i + j == 0:
        return 1
    elif block == 1:
        return 1
    elif block == 0:
        return -1
    else:
        print("Error in decoding matrix!")
        raise

MURAmatrix       = np.zeros([gridSizeX, gridSizeY]);
MURAdecodeMatrix = np.ones([gridSizeX, gridSizeY]);
with open("coded_aperture_array.txt", 'w') as f, open("decoding_matrix.txt", 'w') as d:
    d.write("i,j,d\n") 
    for i in range(0, gridSizeX):
        for j in range(0, gridSizeY):
           
            block  = MURA_code(i, j, gridSizeX);
            decode = MURA_decoding_matrix(i,j, block);
           
            # Write to decoding matrix file
            d.write(str(i) + "," + str(j) + "," + str(decode) + "\n")

            boxLocArrayX = (i) * boxSize - dimX/2. + centeringShift; 
            boxLocArrayY = (j) * boxSize - dimY/2. + centeringShift;

            # Write to Geant geometry construction txt file
            if block == 1:
                f.write(str(round(boxLocArrayX, 4)) + "," 
                        + str(round(boxLocArrayY, 4)) + "\n")
            
                MURAmatrix[i,j] = 1; 

