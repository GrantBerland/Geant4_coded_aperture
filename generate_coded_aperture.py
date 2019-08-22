import numpy as np

boxSize = 0.07;
dimX = 4.;
dimY = 4.;

# p, must be prime!
gridSizeX = 23 
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

def quad_res_mod(i, p):
    '''
    This function checks if i is a quadratic residual of modulo p
    '''
    quad_res_array = [];
    allResidues = np.sqrt(np.mod(range(0, p), p));
    
    for re in allResidues:
        if int(re) == re:
            quad_res_array.append(re);

    if i in quad_res_array:
        return 1;  # is a quadratic residue of mod p
    else:
        return -1; # is a quadratic nonresidue of mod p

def MURA_code(i, j, p):
    '''
    This function generates a Modified Uniform Redundent Aperature
    '''
    if i == 0:
        return 0;
    elif j == 0 and i != 0:
        return 1;
    #elif quad_res_mod(i, p)*quad_res_mod(j, p) == 1:
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

MURAmatrix = np.zeros([gridSizeX, gridSizeY]);
MURAdecodeMatrix = np.ones([gridSizeX, gridSizeY]);
with open("coded_aperture_array.txt", 'w') as f, open("decoding_matrix.txt", 'w') as d:
    d.write("i,j,d\n") 
    for i in range(0, gridSizeX):
        for j in range(0, gridSizeY):
           
            block = MURA_code(i, j, gridSizeX);
            decode = MURA_decoding_matrix(i,j, block);
            
            d.write(str(i) + "," + str(j) + "," + str(decode) + "\n")

            boxLocArrayX = (i) * boxSize; 
            boxLocArrayY = (j) * boxSize;

            if block == 1:
                # Write original and 3 reflections about x, y, and x & y
                f.write(str(round(boxLocArrayX, 4)) + "," 
                        + str(round(boxLocArrayY, 4)) + "\n")
                f.write(str(round(-boxLocArrayX, 4)) + "," 
                        + str(round(boxLocArrayY, 4)) + "\n")
                f.write(str(round(-boxLocArrayX, 4)) + "," 
                        + str(round(-boxLocArrayY, 4)) + "\n")
                f.write(str(round(boxLocArrayX, 4)) + "," 
                        + str(round(-boxLocArrayY, 4)) + "\n")
            
                MURAmatrix[i,j] = 1; 




