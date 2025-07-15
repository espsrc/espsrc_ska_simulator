#!/bin/bash
rm -rf *.MS; ./sorted.sh; python general.py --n_channels=100 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
rm -rf *.MS; ./sorted.sh; python general.py --n_channels=1000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
rm -rf *.MS; ./sorted.sh; python general.py --n_channels=5000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
rm -rf *.MS; ./sorted.sh; python general.py --n_channels=10000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
rm -rf *.MS; ./sorted.sh; python general.py --n_channels=30000 --delta_freq=0.01344 --freq=700 --seconds=600 --cleaning --telescope=SKA-MID-AAstar --catalogue 1 --pixels 1024; ./sorted.sh
