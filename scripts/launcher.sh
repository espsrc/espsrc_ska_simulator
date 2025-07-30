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
python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=3600 --cleaning --telescope=MeerKAT --catalogue 1 --pixels 4096
python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=3600 --cleaning --telescope=SKA1MID --catalogue 1 --pixels 4096
python general.py --n_channels=100 --bandwidth=770 --freq=1283.9 --seconds=3600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 4096
