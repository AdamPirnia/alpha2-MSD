# α₂(t) Pipeline from MD

This repo provides four scripts to extract coordinates from DCD, unwrap PBC jumps,
compute centers of mass, and finally calculate MSD & the non-Gaussian α₂(t).

## Requirements

- Python 3.7+
    + numpy
    + pandas
- VMD 

## Files & Order of Execution

1. **coordinates_extract.py**  
   Extract raw XYZ from DCD using VMD Tcl scripts.

2. **unwrap_coords.py**  
   Remove periodic-boundary wrapping from the raw coordinates.

3. **COM_calc.py**  
   Compute per-frame center-of-mass for each molecule.

4. **alpha2_MSD.py**  
   Compute ⟨|Δr|²⟩, α₂(t) and save to disk.

---

## Usage

All scripts follow the pattern:

## How to run (command example)
```bash
python3 run_pipeline.py \
  --baseDir   /path/to/sim \
  --INdir     trajectories \
  --OUTdir    results \
  --num_dcd   6 \
  --num_mols  1000 \
  --numFrames 1000 \
  --psf       system.psf \
  --dcd       traj.dcd \
  --vmd       /usr/local/bin/vmd \
  --xsc       restart_equil.xsc \
  --num_atoms 3000 \
  --masses    16.00 1.008 1.008