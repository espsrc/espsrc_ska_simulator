# espsrc_ska_simulator
SKA simulations

----- sources.py -------
usage: sources.py [-h] [--path_cfg PATH_CFG] [--path_out PATH_OUT]
                  [--N_srcs N_SRCS] [--rms RMS] [--tofits] [--asksrc]
                  [--backend BACKEND] [--telescope TELESCOPE]

Generate a sky model with sources. If the user does not provide any
parameters, the script will generate a random sky model with 10 sources.

optional arguments:
  -h, --help            show this help message and exit
  --path_cfg PATH_CFG   Path to the config file
  --path_out PATH_OUT   Path to the output file (if this parameter is not
                        provided, then filenames based on the date will be
                        generated)
  --N_srcs N_SRCS       Generate N random sources
  --rms RMS             RMS noise
  --tofits              Save to fits
  --asksrc              Ask for sources
  --backend BACKEND     Imaging backend
  --telescope TELESCOPE
                        Telescope

----- sources.json ------
Example of a sources file
