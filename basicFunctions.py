import os
import matplotlib.image as img 
import matplotlib.pyplot as plt
import numpy as np                                                             
import pandas as pd

def readImage(image):
    """
    Read the .png image and returns a numpy array
    """
    image = img.imread(image)
 
    #This part has to be optimized. It converts a RGB image into a 2D array
    if isinstance(image[0,0], type(np.array([0]))):
        reshaped = np.zeros((len(image), len(image[0])))
        for i in range(len(image)):
            for j in range(len(image[0])):
                reshaped[i,j] = image[i,j,0]
        return reshaped
    else:
        return image


def normalizeImage(imageArray):
    """
    Normalize the image array to minimum value = 0 and maximum value = 1
    """
    
    imageArray = np.divide(imageArray, imageArray.max())
    imageArray[imageArray < 0] = 0
    
    return imageArray

def gradientMagnitude(imageArray, deltaX):
    """
    Calculate the gradient magnitude of a given image array. The pixel length
    or the cell size has to be given to the function. Return the 2D surface
    interface density (Sigma). 
    """
    
    gx, gy = np.gradient(imageArray, deltaX)
    
    return np.hypot(gx, gy)

def curvatureNormalVector(imageArray, deltaX):
    """
    Calculate the mean curvature of a given image array. The curvature (kappa)
    is the negative divergence of the unit vector normal to the surface
    """

    # The normal vector to the surface interface is the gradient of alpha
    gx, gy = np.gradient(imageArray, deltaX)

    # To make the unit vector, the gradient vector is divided by its magnitude
    gradientMag = np.hypot(gx, gy)
    gx /= gradientMag + 0.00000001
    gy /= gradientMag + 0.00000001

    # Then the divergence is calculated
    gxx, gxy = np.gradient(gx, deltaX)
    gyx, gyy = np.gradient(gy, deltaX)

    return -gxx - gyy 

def levelSetFromAlphaField(imageArray, deltaX):
    """
    Calculate the distance function (Level Set) from a pseudo Liquid Volume
    Fraction field (Alpha)
    """
    
    return (2*imageArray-1)*deltaX

def chooseBin(value, serie):
    """
    Assign a bin within a regular serie
    """

    return int((value - serie[0]) // (serie[1] - serie[0])) 

def surfaceCurvatureDistribution(sigma, kappa, alpha, deltaX, 
                                 kMin = -8.5e5, kMax = 8.5e5, nBins = 10000, 
                                 filterSigma = False, filterAlpha = False,
                                 sigmaMin=0.5, alphaMin=0.45, alphaMax=0.55):
    """
    Calculate the Surface Curvature Distribution of a 2pLIF image from sigma,
    kappa and alpha. There is an option to switch on Sigma and Alpha based 
    filters:
        - Sigma Filter: normalized sigma is higher than sigmaMin.
        - Alpha Filter: alpha is between alphaMin and alphaMax

    kMin, kMax are the minimum and maximum curvature values
    nBins is the number of curvature divisions required by the scd. The higher
    the more accurate will be the drop size distribution. 

    Distribution is returned in form of a pandas DataFrame
    """

    print("Calculating the Surface Curvature Distribution")

    #Cretate a new dataFrame that helps to arrange the amount of surface in
    #function of the curvature
    k_bins = [kMin + i * (kMax - kMin) / nBins for i in range(nBins + 1)]
    scd_df = pd.DataFrame([[k, 0] for k in k_bins[:-1]], columns = ["k", "s"])

    #For each pixel the total surface is calculated and arranged inside its 
    #related kappa bin
    for i, (sigma_i, kappa_i, alpha_i) in enumerate(zip(sigma, kappa, alpha)):
        for j, (sigma_ij, kappa_ij, alpha_ij) in enumerate(
        zip(sigma_i, kappa_i, alpha_i)):

            surface = sigma_ij * deltaX * deltaX

            #Check if the data is inside curvature and surface ranges
            insideCurvatureRange = kappa_ij > kMin and kappa_ij < kMax
            insideSurfaceRange = surface > 0 and surface < deltaX
            higherThanSigmaMin = True
            insideAlphaRange = True

            #Apply filters if they are switched on
            if filterSigma: 
                higherThanSigmaMin = sigma_ij * deltaX > sigmaMin
            
            if filterAlpha: 
                higherThanSigmaMin = alpha_ij < alphaMax and alpha_ij > alphaMin
                
            #Surface is then stored into its correspondent bin
            if (insideCurvatureRange and insideSurfaceRange and
                higherThanSigmaMin and insideAlphaRange):
                kappa_bin = chooseBin(kappa_ij, k_bins)
                scd_df.iloc[kappa_bin, scd_df.columns.get_loc("s")] += surface

    return scd_df

def dropSizeDistribution(scd_df, dMin, dMax, nBins, deltaX):
    """
    Calculate the Drop Size Distribution from a Surface Curvature Distribution.
    The SCD must be given in a pandas DataFrame.
    """

    print("Calculating the Drop Size Distribution")

    #The bins and the data frame are created
    d_bins = [dMin + i * (dMax - dMin) / nBins for i in range(nBins + 1)]
    dsd_df = pd.DataFrame([[d, 0] for d in d_bins], columns = ["D", "n"])

    #Calculate the diameter and the number of droplets and place them into
    #the correct diameter bin
    for i, (kappa, surface) in enumerate(zip(scd_df["k"], scd_df["s"])):

        if kappa !=0 and kappa < 4/deltaX:
            diameter = 4/kappa

            if diameter > dMin and diameter < dMax:                                 
                nDrop = surface/np.pi/(4/kappa)**2                                  
                diameter_bin = chooseBin(diameter, d_bins)                                        
                dsd_df.iloc[diameter_bin, dsd_df.columns.get_loc("n")] += nDrop
    
    return dsd_df

def readAllFiles(directory):
    """
    Read all files in a given directory
    """
    return os.listdir(directory)

def makeDir(directory):
    """
    Make a new directory if it does not exist
    """
    if not os.path.exists(directory):                                          
        os.mkdir(directory)    

if __name__ == '__main__':
    image = ""
    imageArray = readImage(image)
    imageArray = normalizeImage(imageArray)
    sigma = gradientMagnitude(imageArray, 2.5e-6)
    kappa = curvatureNormalVector(imageArray, 2.5e-6)
    ls = levelSetFromAlphaField(imageArray, 2.5e-6)
    scd = surfaceCurvatureDistribution(sigma, kappa, imageArray, 2.5e-5, filterSigma = True, filterAlpha = True)
    scd.to_csv("scd.csv")
    dsd = dropSizeDistribution(scd, 0, 80e-6, 40, 2.5e-6)
    dsd.to_csv("dsd.csv")
