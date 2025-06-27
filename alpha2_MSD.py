import numpy as np
import os

######################################################    Parameters    

def a2_MSD(baseDir, INdir, OUTdir, num_dcd, num_mols, numFrames):
    """
    Compute ensemble-averaged mean square displacements (MSD) and non-Gaussian
    parameter α₂(t) for a set of molecules from a series of center-of-mass (COM)
    trajectory files, and save the results.

    Parameters
    ----------
    baseDir : str
        Path to the top-level working directory.
    INdir : str
        Subdirectory under `baseDir` containing the COM data files in
        `com_data/com_<frame>.dat`.
    OUTdir : str
        Subdirectory under `baseDir` where output will be written. Two
        directories will be created inside: `MSDs` and `alpha2s`.
    num_dcd : int
        Number of trajectory (DCD) frames to process. Files are assumed to
        be named `com_0.dat` through `com_{num_dcd-1}.dat`.
    num_mols : int
        Number of molecules whose COM position is stored in each file.
        Each .dat file is expected to have `num_mols*3` columns: x, y, z
        for each molecule.
    numFrames : int
        Minimum number of time-points in each .dat file to include it in the
        average. If a file has fewer lines than `numFrames`, it is skipped
        (and counted as a failure).

    Behavior
    --------
    1. Creates output directories:
         `{baseDir}/{OUTdir}/MSDs`
         `{baseDir}/{OUTdir}/alpha2s`
    2. Loads frame 0 to establish reference positions and accumulates
       ⟨|Δr(t)|²⟩ and ⟨|Δr(t)|⁴⟩ for each subsequent file.
    3. For each frame f in [1, num_dcd):
         - Loads COM data from `{baseDir}/{INdir}/com_data/com_f.dat`.
         - If the file has < numFrames lines, reports and skips it.
         - Otherwise accumulates |r(t)-r(0)|² and its square.
    4. Divides the accumulated sums by (num_dcd - failures) to get the
       ensemble averages.
    5. If num_mols > 2, averages over molecules to produce time-series:
         MSD(t) = ⟨|Δr(t)|²⟩ₜ,ₑₙₛ
         ⟨|Δr(t)|⁴⟩ₜ,ₑₙₛ
    6. Calculates the non-Gaussian parameter:
         α₂(t) = 3⟨|Δr|⁴⟩ / [5 (⟨|Δr|²⟩)²] - 1
    7. Writes out two files:
         `{baseDir}/{OUTdir}/MSDs/MSD_{OUTdir}.dat`
         `{baseDir}/{OUTdir}/alpha2s/a2_{OUTdir}.dat`
       Each file contains one column: the value at each time-step.

    Returns
    -------
    None

    Notes
    -----
    - Make sure that all `.dat` files have the expected shape.
    - `failure` counts how many frames were skipped due to insufficient length.
    - The naming of the output files uses `dir = OUTdir` for consistency.
    """

    os.makedirs(f'{baseDir}/{OUTdir}/MSDs', exist_ok=True)
    os.makedirs(f'{baseDir}/{OUTdir}/alpha2s', exist_ok=True)

    failure = 0

    data = np.loadtxt(f'{baseDir}/{INdir}/com_data/com_0.dat', usecols=range(num_mols*3))
    data = data.reshape(len(data), num_mols, 3)
    data = data - data[0,:,:] 
    data2 = np.einsum('ijk,ijk->ij', data, data)
    data4 = data2**2

    for f in range (1, num_dcd):
        fnam = f'{baseDir}/{dir}/com_data/com_{f}.dat'
        data = np.loadtxt(fnam, usecols=range(num_mols*3))
        data = data.reshape(len(data), num_mols, 3)
        data = data - data[0,:,:] 
        nf = len(data)
        
        if nf < numFrames:
            print (f"len of {fnam} = {len(data)}")
            failure += 1
            continue

        data = np.einsum('ijk,ijk->ij', data, data)
        data2 += data
        data4 += data**2

    data2 /= num_dcd-failure
    data4 /= num_dcd-failure

    if num_mols > 2:                    # Not executed over single water molecule
        data2 = np.mean (data2, axis=1) # (1/num mols) * <|Δr(t)|^2>_t = <|Δr(t)|^2>_t_ens
        data4 = np.mean (data4, axis=1) #                                <|Δr(t)|^4>_t_ens

    alpha2t = 3*data4 / (5*(data2**2)) - 1
    
    np.savetxt(f"{baseDir}/{OUTdir}/MSDs/MSD_{dir}.dat", data2, fmt='%f')
    np.savetxt(f"{baseDir}/{OUTdir}/alpha2s/a2_{dir}.dat", alpha2t, fmt='%f')
