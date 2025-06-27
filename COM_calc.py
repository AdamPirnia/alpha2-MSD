import numpy as np
import pandas as pd
import os

######################################################    Parameters    

def coms(baseDir, INdir, OUTdir, num_dcd, num_mols, num_atoms, particl_mass):
    """
    Compute and save the per-frame center-of-mass coordinates for each molecule
    from unwrapped trajectory snapshots.

    Reads in `num_dcd` files of XYZ coordinates (one file per frame) from
    `{baseDir}/{INdir}/unwrapped/continued_xyz_i.dat`, computes each molecule's
    center of mass, and writes the results to
    `{baseDir}/{OUTdir}/com_data/com_i.dat`.

    Parameters
    ----------
    baseDir : str
        Path to the directory containing both the input and output subdirectories.
    INdir : str
        Name of the subdirectory under `baseDir` where the unwrapped `.dat` files live.
    OUTdir : str
        Name of the subdirectory under `baseDir` where the `com_data` folder will be created.
    num_dcd : int
        Number of trajectory frames (i.e., number of input files to process).
    num_mols : int
        Number of molecules per frame.
    num_atoms : int
        Number of atoms in each molecule.
    particl_mass : list of float
        Atomic masses for the `num_atoms` atoms in a single molecule,
        in the same order as the coordinates appear in the input files.
    stride : int, optional
        Frame-skipping stride (not yet implemented; defaults to 1).

    Returns
    -------
    None
        This function writes one output file per frame and does not return a value.

    Side effects
    ------------
    - Suppresses FutureWarnings.
    - Creates directory `{baseDir}/{OUTdir}/com_data` if it does not exist.
    - Prints each filename as it is processed.
    - Saves flattened center-of-mass arrays to text files.

    Notes
    -----
    - Input files must be whitespace-delimited XYZ snapshots, with each line
      giving one atom's x,y,z in sequence.  The total number of lines per file
      must equal `num_mols * num_atoms`.
    - Output files are whitespace-delimited, with each row corresponding to one
      molecule's (x, y, z) center of mass.

    Example
    -------
    >>> coms(
    ...     baseDir="/home/user/sim",
    ...     INdir="traj",
    ...     OUTdir="analysis",
    ...     num_dcd=1000,
    ...     num_mols=500,
    ...     num_atoms=3,
    ...     particl_mass=[16.00, 1.008, 1.008]  # e.g. water: O, H, H
    ... )
    """
    warnings.simplefilter(action='ignore', category=FutureWarning)

    def cmass(n, frames, masses, mt):
        """
        Compute centers of mass for each molecule in each frame.

        Parameters
        ----------
        n : int
            Number of molecules.
        frames : ndarray, shape (n_frames, n, num_atoms*3)
            Atom coordinates flattened per molecule.
        masses : ndarray, shape (num_atoms,)
            Atomic masses for one molecule.
        mt : float
            Total mass of one molecule (sum of `masses`).

        Returns
        -------
        com : ndarray, shape (n_frames, n, 3)
            Center-of-mass coordinates for each molecule in each frame.
        """
        # Reshape the input data into (frames, molecules, atoms, 3) for easier computation
        frames = np.reshape(frames, (frames.shape[0], n, num_atoms, 3))
        # Calculate the center of mass for each water molecule in every frame
        com = np.sum(frames * masses[:, np.newaxis], axis=2) / mt
        return com
    
    os.makedirs(f'{baseDir}/{OUTdir}/com_data', exist_ok=True)   
    for i in range (num_dcd):
        total_mass = np.sum(np.array(particl_mass))
        fnam = f'{baseDir}/{INdir}/unwrapped/continued_xyz_{i}.dat'    
        print (f"{fnam}")
        data = pd.read_csv(fnam, delim_whitespace=True, header=None).values

        # np.savetxt(f'{dir}/results/var_{i}ns.dat', vari))
        data = np.reshape(data, (np.shape(data)[0], num_mols, num_atoms*3)) 

        centers_of_mass = cmass(num_mols, data, particl_mass, total_mass)

        # Reshape the centers_of_mass array to (1000, num_mol * 3)
        centers_of_mass_flattened = centers_of_mass.reshape(data.shape[0], -1) 
        np.savetxt(f'{baseDir}/{OUTdir}/com_data/com_{i}.dat', centers_of_mass_flattened)
