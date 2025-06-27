import numpy as np
import os

    
def unwrapper (baseDir, INdir, OUTdir, xsc, num_dcd, num_atoms, stride=1):
    os.makedirs(f'{baseDir}/{OUTdir}/unwrapped', exist_ok=True)
    """
    Read a series of coordinates directly extractyed from MD trajectories, unwrap 
    periodic boundary crossings, and write out unwrapped coordinate files.

    This function:
      1. Reads the last line of the `.xsc` file to determine the simulation box vectors.
      2. For each frame in the input directory:
         a. Loads the wrapped coordinates from `xyz_i.dat`.
         b. Computes frame-to-frame displacements.
         c. Applies minimum-image corrections (i.e. if a jump > half box, subtract/add
            full box length).
         d. Cumulatively sums the corrected displacements to produce unwrapped
            trajectories.
         e. Saves the unwrapped coordinates to the output directory.

    Parameters
    ----------
    baseDir : str
        Root directory for the simulation analysis. The XSC, INdir, and OUTdir
        are all relative to this.
    INdir : str
        Subdirectory under `baseDir` containing the wrapped coordinate files
        named `xyz_0.dat, xyz_1.dat, …, xyz_{num_dcd-1}.dat`.
    OUTdir : str
        Subdirectory under `baseDir` where unwrapped files will be written.
        Files will be named `unwrapped_xyz_0.dat`, etc. Created if needed.
    xsc : str
        Path (relative to `baseDir`) to the `.xsc` file holding box vectors.
        The last line is expected to contain at least 9 numeric entries
        from which the x-, y-, and z-box lengths are extracted.
    num_dcd : int
        Number of frames to process (i.e. number of `xyz_i.dat` files).
    num_atoms : int
        Number of atoms per frame. Each `xyz_i.dat` is expected to have
        `num_atoms*3` columns when flattened.
    stride : int, optional
        Frame stride for processing. Currently reserved (not implemented);
        defaults to 1 (i.e. process every frame).

    Notes
    -----
    - The box lengths are taken from indices [1], [5], [9] of the last XSC line.
    - Half-box thresholds (`thold`) are used for minimum-image unwrapping.
    - The function reshapes each flat XYZ array into shape (n_frames, num_atoms, 3)
      before unwrapping, and flattens it back when saving.
    - Outputs are written in plain text with six-decimal floats.

    Examples
    --------
    >>> unwrapper(
    ...     baseDir="/home/user/sim",
    ...     INdir="anlz/NVT/wrapped",
    ...     OUTdir="anlz/NVT",
    ...     xsc="restart_equil.xsc",
    ...     num_dcd=1000,
    ...     num_atoms=1234
    ... )
    """

    with open (f"{baseDir}/{xsc}", 'r') as fr:
        lins = fr.readlines()
    box_size = np.array ([lins[-1].split()[1], lins[-1].split()[5], lins[-1].split()[9]]).astype(float)
    thold = box_size/2

    for i in range(num_dcd):
        anlz_file = f'{baseDir}/{INdir}/xyz_{i}.dat'
        xyz = np.loadtxt(anlz_file)[::stride]
        xyz = xyz.reshape(len(xyz), num_atoms, 3)
        diffs = xyz[1:] - xyz[:-1]
        new_xyz = np.zeros_like(xyz)
        new_xyz[0] = xyz[0]
        
        # Correct differences along each dimension
        for dim in range(3):
            diffs[..., dim][diffs[..., dim] > thold[dim]] -= box_size[dim]
            diffs[..., dim][diffs[..., dim] < -thold[dim]] += box_size[dim]
        
        # Compute cumulative sum once after all corrections
        new_xyz[1:] = new_xyz[0] + np.cumsum(diffs, axis=0)
        
        shp = np.shape(new_xyz)
        print(shp)
        np.savetxt(f"{baseDir}/{OUTdir}/unwrapped_xyz_{i}.dat", new_xyz.reshape(shp[0], num_atoms*3), fmt='%f')
        
    del xyz
    del new_xyz
