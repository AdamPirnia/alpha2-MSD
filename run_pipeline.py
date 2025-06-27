#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os

def call_script(script, args):
    """Helper to call one of your scripts with the given args."""
    cmd = [sys.executable, os.path.join(os.getcwd(), script)] + args
    print(f"\n→ Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    p = argparse.ArgumentParser(
        description="One-shot pipeline: extract → unwrap → COM → alpha2/MSD"
    )
    # common args
    p.add_argument("--baseDir",   required=True, help="Top-level directory")
    p.add_argument("--INdir",     required=True, help="Input subdir under baseDir")
    p.add_argument("--OUTdir",    required=True, help="Output subdir under baseDir")
    p.add_argument("--num_dcd",   type=int, required=True, help="Number of DCD frames")
    p.add_argument("--num_mols",  type=int, required=True, help="Number of molecules")
    p.add_argument("--numFrames", type=int, required=True,
                   help="Min frames per file for MSD step")

    # step-specific args
    p.add_argument("--psf",       required=True, help="PSF filename (for VMD step)")
    p.add_argument("--dcd",       required=True, help="DCD filename (for VMD step)")
    p.add_argument("--vmd",       required=True, help="Path to VMD executable")
    p.add_argument("--xsc",       required=True, help="XSC restart file (unwrap step)")
    p.add_argument("--num_atoms", type=int, required=True,
                   help="Number of atoms (unwrap/COM steps)")
    p.add_argument("--masses",    nargs="+", type=float, required=True,
                   help="List of atomic masses (COM step)")

    args = p.parse_args()

    common = [
        "--baseDir",   args.baseDir,
        "--INdir",     args.INdir,
        "--OUTdir",    args.OUTdir,
        "--num_dcd",   str(args.num_dcd),
        "--num_mols",  str(args.num_mols),
    ]

    # 1) coordinates_extract.py
    step1 = common + [
        "--psf", args.psf,
        "--dcd", args.dcd,
        "--vmd", args.vmd,
    ]
    call_script("coordinates_extract.py", step1)

    # 2) unwrap_coords.py
    step2 = common + [
        "--xsc",       args.xsc,
        "--num_atoms", str(args.num_atoms),
    ]
    call_script("unwrap_coords.py", step2)

    # 3) COM_calc.py
    step3 = common + [
        "--num_atoms", str(args.num_atoms),
        "--masses",   *map(str, args.masses),
    ]
    call_script("COM_calc.py", step3)

    # 4) alpha2_MSD.py
    step4 = common + [
        "--numFrames", str(args.numFrames),
    ]
    call_script("alpha2_MSD.py", step4)

if __name__ == "__main__":
    main()
