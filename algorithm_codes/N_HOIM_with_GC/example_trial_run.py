#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

# ---------------------- TWEAK THESE ----------------------
TARGET_SCRIPT = "run_one_trial.py"   
PREFIX = "uf20"
INSTANCE = 2
THLD_DELTA = "0.001"
MAX_ITER = 50000          # e.g., 5000, or keep None to use default formula inside the script
TRIAL_IDS = [1]    # sample runs
# ---------------------------------------------------------

def main():
    base_dir = Path(__file__).resolve().parent
    script_path = base_dir / TARGET_SCRIPT
    if not script_path.exists():
        raise FileNotFoundError(f"Cannot find {script_path}. Update TARGET_SCRIPT.")

    for tid in TRIAL_IDS:
        cmd = [
            sys.executable, str(script_path),
            "--prefix", str(PREFIX),
            "--instance", str(INSTANCE),
            "--thld-delta", str(THLD_DELTA),
            "--trial-id", str(tid),
        ]
        if MAX_ITER is not None:
            cmd += ["--max-iter", str(MAX_ITER)]

        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
