import subprocess


path = "/mnt/scratch/espsrc_ska_simulator/scripts/20250724_1216/20250724_1216_visibilities.MS"
command = f"OPENBLAS_NUM_THREADS=1 wsclean -weight briggs 0.0 -multiscale -size 2048 2048 -scale 0.0009745046474780809deg -niter 5000 -mgain 0.8 -auto-threshold 0.3 -auto-mask 3 -channels-out 8 -join-channels ${path}"

completed_process = subprocess.run(
    command,
    shell=True,
    capture_output=True,
    text=True,
    # Raises exception on return code != 0
    check=True,
)