import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import csv

from scipy import signal
from scipy import ndimage
from scipy import misc

import cv2 


def makeMURA(gridSizeX=17):
    detSize = 4.


    # Creates MURA coded aperture and decoding matrices
    boxSize = 0.22;
    dimX = 4.;
    dimY = 4.;

    centeringShift = 0.22;

    # p, must be prime and satify L = 4m + 1
    gridSizeY = gridSizeX; # (square)

    def jacobi_symbol(a, n):
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
        #
        #This function generates a Modified Uniform Redundent Aperature
        #
        if i == 0:
            return 0;
        elif j == 0 and i != 0:
            return 1;
        elif jacobi_symbol(i, p)*jacobi_symbol(j, p) == 1: 
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
    MURAdecodeMatrix = -np.ones([gridSizeX, gridSizeY]);
    with open("coded_aperture_array.txt", 'w') as f, open("decoding_matrix.txt", 'w') as d:
        d.write("i,j,d\n") 
        for i in range(0, gridSizeX):
            for j in range(0, gridSizeY):

                block = MURA_code(i, j, gridSizeX);
                decode = MURA_decoding_matrix(i,j, block);
                MURAdecodeMatrix[i,j] = decode;

                d.write(str(i) + "," + str(j) + "," + str(decode) + "\n")

                boxLocArrayX = (i) * boxSize + centeringShift; 
                boxLocArrayY = (j) * boxSize + centeringShift;

                if block == 1:
                    f.write(str(round(boxLocArrayX, 4)) + "," 
                            + str(round(boxLocArrayY, 4)) + "\n")

                    MURAmatrix[i,j] = 1;
                    
    return MURAmatrix, MURAdecodeMatrix


def getMeanImage(hits, totalDecode, MURAdecodeMatrix, plotOn, mask_info):
    
    def retrieveSignal(det, decoder):
        sig = np.fliplr(signal.fftconvolve(det, 
                             decoder, 
                             mode='full'))
        return sig
    
    nBins   = 16
    limits  = 4.1

    pixelSize    = 2.2/10;   # cm
    pixelSpacing = 0.26/10;  # cm

    pixelCoordsX = [];
    pixelCoordsY = [];

    detectorSpacing = 20.5/10; # cm
    centering       = 0.32/10  # cm

    detectorOffsetX = [detectorSpacing, -detectorSpacing, detectorSpacing, -detectorSpacing];
    detectorOffsetY = [detectorSpacing, -detectorSpacing, -detectorSpacing, detectorSpacing];

    for detectorNum in range(0, 4):
        for i in range(0, nBins+1):
            pixelCoordsX.append((pixelSize+pixelSpacing)*(i-nBins/2) + detectorOffsetX[detectorNum] + centering);
            pixelCoordsY.append((pixelSize+pixelSpacing)*(i-nBins/2) + detectorOffsetY[detectorNum] + centering);


    pixelCoordsX = np.unique(np.sort(pixelCoordsX))
    pixelCoordsY = np.unique(np.sort(pixelCoordsY))

    h_sig, xedges, yedges, im = plt.hist2d(hits.x, hits.y, bins=(pixelCoordsX, pixelCoordsY));

    deconvolved_sig = signal.fftconvolve(h_sig,
                                         totalDecode, 
                                         mode='full');

    det1 = h_sig[0:nBins  , 0:nBins]    # top left
    det2 = h_sig[0:nBins  , nBins:-1]   # top right
    det3 = h_sig[nBins:-1 , 0:nBins]    # bottom left
    det4 = h_sig[nBins:-1 , nBins:-1]   # bottom right


    det1 -= np.mean(det1)
    det2 -= np.mean(det2)
    det3 -= np.mean(det3)
    det4 -= np.mean(det4)

    sig1 = retrieveSignal(det1, MURAdecodeMatrix)
    sig2 = retrieveSignal(det2, MURAdecodeMatrix)
    sig3 = retrieveSignal(det3, MURAdecodeMatrix)
    sig4 = retrieveSignal(det4, MURAdecodeMatrix)

    sig1[sig1 < 0] = 0
    sig2[sig2 < 0] = 0
    sig3[sig3 < 0] = 0
    sig4[sig4 < 0] = 0


    mask       = np.zeros([sig1.shape[0], sig1.shape[1]])
    mean_image = np.zeros([sig1.shape[0], sig1.shape[1]])
    mean_image += sig1
    mean_image += sig2
    mean_image += sig3
    mean_image += sig4
    mean_image /= 4;

    dO = 720;
    A   = mask_info[0];
    shX = mask_info[1];
    shY = mask_info[2];

    for i in range(0, dO):
        mask[(A*np.cos(dO*i) + shX).astype(np.int64), 
                (A*np.sin(dO*i) + shY).astype(np.int64)] = 1;
        for j in range(1, A+1):
            mask[((A-j)*np.cos(dO*i) + shX).astype(np.int64),
                ((A-j)*np.sin(dO*i) + shY).astype(np.int64)] = 2;

    if plotOn == 1:
        plt.figure(figsize=(9,4))
        plt.subplot(1,2,1);
        plt.colorbar(im);
        lims = 4.2
        plt.xlim([-lims , lims])
        plt.ylim([-lims , lims])
    
        plt.subplot(1,2,2);
        mapable = plt.imshow(mean_image, origin='left')
        plt.colorbar(mapable);
    else:
        plt.clf();
    return mean_image, mask

def plot3DmeanImages(mean_image):
    smoothing_factor = 5
    kernel = np.ones((smoothing_factor,smoothing_factor),np.float32)/(smoothing_factor**2)
    dst = cv2.filter2D(mean_image,-1,kernel)


    blur_factor = 5
    dst_blur = cv2.GaussianBlur(mean_image, (blur_factor,blur_factor),0)


    x = np.linspace(0, len(mean_image), len(mean_image))
    y = np.linspace(0, len(mean_image), len(mean_image))

    X, Y = np.meshgrid(x, y)

    fig = plt.figure(3, figsize=(12,16))
    ax = fig.add_subplot(321, projection='3d')
    ax.plot_surface(X, Y, 
                    mean_image, 
                    cmap='binary');
    plt.title('Original Deconvolution')
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.subplot(322)
    #plt.imshow(mean_image, origin='left', extent=[0,32,0,32]);
    plt.pcolormesh(mean_image)
    plt.title('Original Deconvolution')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar();

    ax = fig.add_subplot(323, projection='3d')
    ax.plot_surface(X, Y, 
                    dst, 
                    cmap='binary', 
                    alpha=0.5);
    plt.title('Smoothed Deconvolution');
    plt.xlabel('X')
    plt.ylabel('Y')

    max_dst = np.unravel_index(np.argmax(dst), dst.shape)
    max_dst_val = np.amax(dst)


    plt.plot(max_dst[0]*np.ones(10), max_dst[1]*np.ones(10), np.linspace(-30, max_dst_val, 10), 
              markersize=5, color='r')

    plt.subplot(324)
    #plt.imshow(dst, origin='left', extent=[0,32,0,32]);
    plt.pcolormesh(dst)
    plt.title('Smoothed Deconvolution');
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar();

    dst_max = np.unravel_index(np.argmax(dst), dst.shape)
    plt.gca().add_patch(patches.Circle((dst_max[0], dst_max[1]), radius=1., 
                                       facecolor=None,
                                       linewidth=1,
                                       edgecolor='r',
                                       fill=False));


    ax = fig.add_subplot(325, projection='3d')
    ax.plot_surface(X, Y, 
                    dst_blur, 
                    cmap='binary');
    plt.title('Gaussian Blurred Deconvolution');
    plt.xlabel('X')
    plt.ylabel('Y')

    plt.subplot(326)
    #plt.imshow(dst_blur, origin='left', extent=[0,32,0,32]);
    plt.pcolormesh(dst_blur);
    plt.title('Gaussian Blurred Deconvolution');
    plt.xlabel('X')
    plt.ylabel('Y');
    plt.colorbar();

    dst_blur_max = np.unravel_index(np.argmax(dst_blur), dst_blur.shape)
    plt.gca().add_patch(patches.Circle((dst_blur_max[0], dst_blur_max[1]), radius=1., 
                                       facecolor=None,
                                       linewidth=1,
                                       edgecolor='r',
                                       fill=False));

def plotAllMeanImages(names, *mean_images):
    
    N = len(mean_images[0]);
    
    for i in range(0, N):
        
        mean_image = mean_images[0][i]
        
        blur_factor = 5
        dst_blur = cv2.GaussianBlur(mean_image, (blur_factor,blur_factor),0)

        x = np.linspace(0, len(mean_image), len(mean_image))
        y = np.linspace(0, len(mean_image), len(mean_image))

        X, Y = np.meshgrid(x, y)

        fig = plt.figure(i, figsize=(12,12))
        plt.suptitle(names[i])
        
        ax = fig.add_subplot(221, projection='3d')
        ax.plot_surface(X, Y, 
                        mean_image, 
                        cmap='binary');
        plt.title('Original Deconvolution')
        plt.xlabel('X')
        plt.ylabel('Y')

        plt.subplot(222)
        plt.pcolormesh(mean_image)
        plt.title('Original Deconvolution')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.colorbar();

        ax = fig.add_subplot(223, projection='3d')
        ax.plot_surface(X, Y, 
                        dst_blur, 
                        cmap='binary');
        plt.title('Gaussian Blurred Deconvolution');
        plt.xlabel('X')
        plt.ylabel('Y')

        plt.subplot(224)
        plt.pcolormesh(dst_blur);
        plt.title('Gaussian Blurred Deconvolution');
        plt.xlabel('X')
        plt.ylabel('Y');
        plt.colorbar();

        dst_blur_max = np.unravel_index(np.argmax(dst_blur), dst_blur.shape)
        plt.gca().add_patch(patches.Circle((dst_blur_max[0], dst_blur_max[1]), radius=1., 
                                           facecolor=None,
                                           linewidth=1,
                                           edgecolor='r',
                                           fill=False));  
   

def twoD_Gaussian(x, y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel() 


