#!/bin/bash
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=100 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=1000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=5000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=10000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=30000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh

# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=100 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=1000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=5000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=10000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 1024; ./sorted.sh
# rm -rf *.MS; ./sorted.sh; python general.py --n_channels=30000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 1024; ./sorted.sh

# for i in  1000 5000 10000 30000; do
#     ./sorted.sh
#     python general.py --n_channels=${i} --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096 --fov 1
#     ./sorted.sh
# done

# 20250729
#python general.py --n_channels=100 --bandwidth=100 --freq=700 --seconds=600 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 2048 
#python general.py --n_channels=100 --bandwidth=100 --freq=700 --seconds=600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 2048
#python general.py --n_channels=100 --bandwidth=100 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 2048 

# 20250730
#python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=3600 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
#python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=3600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
#python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=3600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096

#20250801
# python general.py --n_channels=1000 --bandwidth=770 --freq=1283.9 --seconds=36000 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
# python general.py --n_channels=1000 --bandwidth=770 --freq=1283.9 --seconds=36000 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
# python general.py --n_channels=1000 --bandwidth=770 --freq=1283.9 --seconds=36000 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096
#20250816
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=36000 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=36000 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=36000 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096

#2025-09-8
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=819.2 --freq=1265.6 --seconds=8 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=819.2 --freq=1265.6 --seconds=8 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=819.2 --freq=1265.6 --seconds=8 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096
# 
#2025-09-09
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=819.2 --freq=1265.6 --seconds=8 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=720 --freq=1310 --seconds=8 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=720 --freq=1310 --seconds=8 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096

# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=819.2 --freq=1265.6 --seconds=80 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=720 --freq=1310 --seconds=80 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=720 --freq=1310 --seconds=80 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096
# 
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=819.2 --freq=1265.6 --seconds=800 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=720 --freq=1310 --seconds=800 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
# find . -name *.MS -exec rm -rf {} \;
# python general.py --n_channels=4096 --bandwidth=720 --freq=1310 --seconds=800 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096


# Band 1
#python general.py --catalogue=1 --fov=90arcsec --pixels=4096 --seconds=600 --cleaning --bandwidth=700 --freq=700 --n_channels=2170 --telescope=SKA1MID --center="10h00m27.4474s,+02d20m57s"

 python general.py --fits=./MIGHTEE_Continuum_DR1_COSMOS_5p2arcsec_I_v1.1_FinalCatalogue.srl.fits --pixels=1024 --seconds=600 --cleaning --bandwidth=700 --freq=700 --n_channels=2170 --telescope=SKA1MID --center="10h00m27.4474s,+02d20m57s" --column-mapping="0,2,4,6,-1,-1,-1,-1,-1,-1,10,12,14"
