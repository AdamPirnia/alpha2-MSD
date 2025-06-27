import os
import subprocess

def raw_coords (baseDir, INdir, OUTdir, psf, dcd, num_dcd, num_mols, vmd):
    """
    Extracts raw XYZ coordinates from a series of DCD trajectories using VMD,
    by generating per-segment Tcl scripts, running them in batch, and saving
    the results as text files.

    For each trajectory chunk (0→1 ns, 1→2 ns, …), this function:
      1. Creates output directories (`OUTdir/data`, `writenCodes`, `logs`).
      2. Writes a VMD Tcl script that:
         - Loads the PSF and DCD for that chunk.
         - Iterates over every frame.
         - Selects residues 0 through `num_mols` (e.g. water molecules).
         - Extracts their x,y,z coordinates.
         - Writes one line per frame to `OUTdir/data/xyz_{i}.dat`.
      3. Invokes VMD in text mode (`-dispdev text -nt 64`) to run the script.
      4. Captures stdout/stderr to `logs/log_{i}.lgo`.

    Parameters
    ----------
    baseDir : str
        Root directory containing the `INdir` subfolders named
        “0to1ns”, “1to2ns”, … up to `num_dcd-1`to`num_dcd`ns.
    INdir : str
        Name of the subfolder under `baseDir` where all trajectory chunks live.
    OUTdir : str
        Name of the output subdirectory under `baseDir` in which `data/` will be 
        created to hold the `.dat` files. OUTdir should exist.
    psf : str
        Base filename (without extension) of the topology file in each chunk
        (e.g. `"system"` if your files are `system.psf`).
    dcd : str
        Base filename (without extension) of the trajectory file in each chunk
        (e.g. `"traj"` if your files are `traj.dcd`).
    num_dcd : int
        Number of trajectory chunks to process.  Generates scripts for
        indices `0` through `num_dcd-1`.
    num_mols : int
        Maximum residue index to select (e.g. number of water molecules);
        selects residues `0` through `num_mols`.
    vmd :
        The directory of the executable of the VMD software.
    
    Side Effects
    ------------
    - Creates directories:
        - `{baseDir}/{OUTdir}/data`
        - `writenCodes`
        - `logs`
    - Writes Tcl scripts to `writenCodes/coords_{i}.tcl`.
    - Runs VMD on each script, logging to `logs/log_{i}.lgo`.
    - Produces coordinate files `xyz_{i}.dat` in `{baseDir}/{OUTdir}/data`.

    Returns
    -------
    None

    Example
    -------
    >>> raw_coords(
    ...     baseDir="/home/user/sim",
    ...     INdir="trajectories",
    ...     OUTdir="extracted",
    ...     psf="system",
    ...     dcd="traj",
    ...     num_dcd=6,
    ...     num_mols=1000
    ...     /scratch/user/programs/vmddir/bin/vmd
    ... )
    # Creates:
    #   /home/user/sim/extracted/data/xyz_0.dat, …, xyz_5.dat
    #   writenCodes/coords_0.tcl, …, coords_5.tcl
    #   logs/log_0.lgo, …, log_5.lgo
    """

    os.makedirs(f'{baseDir}/{OUTdir}/data', exist_ok=True)
    os.makedirs(f'writenCodes', exist_ok=True)
    os.makedirs(f'logs', exist_ok=True)

    for i in range (num_dcd):
        with open (f"writenCodes/coords_{i}.tcl", "w") as f:
            f.write(f"""
# Print the start message with date and time
puts "Starting Analysis - date/time"                                            
date                                                                            


# Set the base directory path to the trajectories
set baseDir {baseDir}
set resDir {baseDir}/{OUTdir}

set dir {INdir}

set outfile [open  ${{resDir}}/data/xyz_{i}.dat w]
set PSF "${{baseDir}}/${{dir}}/{i}to{i+1}ns/{psf}.psf"
set DCD "${{baseDir}}/${{dir}}/{i}to{i+1}ns/{dcd}.dcd"

# Load the velocity trajectory file (PSF and DCD files)
set traj [mol load psf $PSF dcd $DCD]
puts "Trajectory is loaded"

# Get the number of frames in the trajectory
set nf [molinfo $traj get numframes]  
puts $nf 

# Loop over each frame in the trajectory
for {{set frame 0}} {{$frame < $nf}} {{incr frame}} {{
    # Go to the current frame
    animate goto $frame

    # Select water molecules
    set sel [atomselect top "residue 0 to {num_mols}"]

    set coor [$sel get {{x y z}}]

    set frameData ""
    # Loop over each set of 3 coordinates in the coor list
    for {{set j 0}} {{$j < [llength $coor]}} {{incr j}} {{
        
        set x [expr {{double([lindex [lindex $coor $j] 0])}}]
        set y [expr {{double([lindex [lindex $coor $j] 1])}}]
        set z [expr {{double([lindex [lindex $coor $j] 2])}}]


        # Print the coordinates with 3 decimal places into the output file
        append frameData [format " %.6f %.6f %.6f " $x $y $z]
    }}
    puts $outfile $frameData
}}

mol delete all
close $outfile

#garbage collection
gc
exit
""")
        
        print ("after writing file.")
        # Run the bash command using subprocess
        command1 = f"{vmd} -dispdev text -nt 64 -e writenCodes/coords_{i}.tcl > logs/log_{i}.lgo"
        print ("after executing the command.")
        process = subprocess.run(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stdout = process.stdout.decode('utf-8')
        stderr = process.stderr.decode('utf-8')

        print(f"STDOUT:\n{stdout}")
        print(f"STDERR:\n{stderr}")
    
