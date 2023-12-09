import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.io
plt.jet()

# Vonk 3d function from SGems
def vonk3d(rseed, dx, dy, dz, ax, ay, az, ix, iy, iz, pop, med, nu, vel=[1]):
    """
    Generates a 3D random field based using von Karman Dist.
    
    Parameters:
    -----------
    rseed : int
        Random seed for reproducibility.
    dx, dy, dz : float
        Grid spacing in the x, y, and z directions, respectively.
    ax, ay, az : float
        Anisotropy parameters in the x, y, and z directions, respectively.
    ix, iy, iz : int
        Size of the domain in the x, y, and z directions, respectively.
    pop : int
        Population type. Can be either:
        1 - Gaussian
        2 - PDF
    med : int
        Medium type determining the correlation function:
        1 - Gaussian
        2 - Exponential
        3 - von Karman
    nu : float
        Parameter related to the von Karman correlation function.
    vel : list, optional
        List containing velocity values to scale the generated random field. Default is [1].
        
    Returns:
    --------
    ndarray
        A 3D array representing a normalized random field with zero mean and a standard deviation of 1.
    """

    if pop == 1:  # POP=GAUSSIAN
        D = 3 - nu
    elif pop == 2:  # POP=PDF
        D = 3 - nu / 2

    nx = np.ceil(ix / dx)
    ny = np.ceil(iy / dy)
    nz = np.ceil(iz / dz)

    # SPACEDOMAIN GRID
    x = np.linspace(1, int(nx), int(nx)) * dx
    y = np.linspace(1, int(ny), int(ny)) * dy
    z = np.linspace(1, int(nz), int(nz)) * dz

    # WAVENUMBER DOMAIN GRID
    dkx = 1 / (nx * dx)
    dky = 1 / (ny * dy)
    dkz = 1 / (nz * dz)
    [kx, ky, kz] = np.meshgrid(2 * math.pi * dkx * np.linspace(int(-nx / 2),
                                                               int(nx / 2 - 1),
                                                               int(nx / 2 + nx / 2 - 1 + 1)),
                               2 * math.pi * dky * np.linspace(int(-ny / 2),
                                                               int(ny / 2 - 1),
                                                               int(ny / 2 + ny / 2 - 1 + 1)),
                               2 * math.pi * dkz * np.linspace(int(-nz / 2),
                                                               int(nz / 2 - 1),
                                                               int(nz / 2 + nz / 2 - 1 + 1)))

    k = np.sqrt(kx ** 2 * ax ** 2 + ky ** 2 * ay ** 2 + kz ** 2 * az ** 2)

    if med == 1:
        # Gaussian Chh
        # SQRT of the FOURIER TRANSROMF OF Gaussian CORR fkt.
        expcorr = ((ax * az * ay) / 2) * np.exp(-(k ** 2) / 2)  # (Exact F(C_hh(x))
        expcorr = expcorr / np.amax(expcorr)  # normalizing sqr(F(C_hh))
        expcorr = np.sqrt(expcorr)

    if med == 2:
        # Exponential Chh
        # SQRT of the FOURIER TRANSROMF OF exp CORR fkt.
        expcorr = 1 / ((1 + k ** 2) ** (1.5))  # (Exact F(C_hh(x))
        expcorr = expcorr / np.amax(expcorr)  # normalizing sqr(F(C_hh))
        expcorr = np.sqrt(expcorr);

    if med == 3:
        # von Karman
        # SQRT of the FOURIER TRANSROMF OF vonk CORR fkt.
        expcorr = 1 / ((1 + k ** 2) ** (nu + 1));  # (Exact F(C_hh(x))
        expcorr = expcorr / np.amax(expcorr);  # normalizing sqr(F(C_hh))
        expcorr = np.sqrt(expcorr);  #

    # DATA
    #     rng = np.random.default_rng(seed=rseed)
    np.random.seed(rseed)
    data = np.random.random(size=(int(ny), int(nx), int(nz)))

    # GOING TO FOURIER DOMAIN
    fdata = np.fft.fftshift(np.fft.fftn(data))

    # MULTIPLYING fdata by sqrt(C_hh)
    newfdata = fdata * expcorr

    # FROM FOURIER TO SPACE DOMAIN
    randdata = np.real(np.fft.ifftn(np.fft.fftshift(newfdata)))
    rdata = randdata;

    if pop == 1:
        # scaling filed according to vel
        rdata = randdata * 2 * vel[0]
        data = data * 2 * vel[0]

    return (rdata - rdata.mean()) / rdata.std()


def scale_random_field(kmap, mean, std):
    # kmap is a normal distribution with (0, 1)
    # this function remaps such a normal distribution to a lognormal distribution with specified mean of logK 
    logk_mean = np.log( mean**2 / np.sqrt(mean**2 + std**2))
    print(logk_mean)
    logk_var = np.log( 1 + std**2 / mean**2)
    print(logk_var)
    rescaled_logk = kmap * np.sqrt(logk_var) + logk_mean
    
    return np.exp(rescaled_logk)

nX = 30;
nY = 30;
nZ = 260;
mean = 75e-15 # 75 mD
std = 25e-13

ax = np.random.randint(1*nX//2, 3*nX//4, (1, 1))
ay = np.random.randint(1*nY//2, 3*nY//4, (1, 1))
az = np.random.randint(1*nZ//2, 3*nZ//4, (1, 1))
kmap = np.zeros((nX, nY, nZ))
kmap = vonk3d(rseed=1,
              dx=1, dy=1, dz=1, 
              ax=ax[0][0], ay=ay[0][0], az=az[0][0],
              ix=nX, iy=nY, iz=nZ,
              pop=1, med=3,
              nu=1, vel=[1])

rescaled_kmap = scale_random_field(kmap, mean, std);

minK = 1e-15;
maxK = 1e-12;
for k in range(0,nZ):
    perm = rescaled_kmap[:,:,k]
    logPerm = np.log(perm)
    logPerm = (logPerm-np.min(logPerm))/(np.max(logPerm)-np.min(logPerm))*(np.log(maxK)-np.log(minK))+np.log(minK)
    perm = np.exp(logPerm)
    rescaled_kmap[:,:,k] = perm

plt.figure()
plt.imshow(np.log10(rescaled_kmap[:,:,0]))
plt.colorbar()
plt.axis('equal')
plt.savefig('test0.pdf')
    
rescaled_kmap = np.reshape(rescaled_kmap, (nX*nY*nZ,1),order='F')
np.savetxt('permeability30x30x1000.out', rescaled_kmap, delimiter=' ')


# plt.figure()
# plt.imshow(np.log10(rescaled_kmap[:,:,1]))
# plt.colorbar()
# plt.axis('equal')
# plt.savefig('test1.pdf')

# plt.figure()
# plt.imshow(np.log10(rescaled_kmap[:,:,2]))
# plt.colorbar()
# plt.axis('equal')
# plt.savefig('test2.pdf')

# plt.figure()
# plt.imshow(np.log10(rescaled_kmap[:,:,3]))
# plt.colorbar()
# plt.axis('equal')
# plt.savefig('test3.pdf')

# plt.figure()
# plt.imshow(np.log10(rescaled_kmap[:,:,4]))
# plt.colorbar()
# plt.axis('equal')
# plt.savefig('test4.pdf')
