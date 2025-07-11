import matplotlib.pyplot as plt
import numpy as np
import csv
import argparse
from utils import show_exc
import datetime, os, sys

def plot_monitor_data(input_file, output_file=None, rowheader=False):
    """
    Plots data from a CSV file and saves it as an image.
    :param input_file: Path to the input CSV file.
    :param output_file: Path to the output image file.
    """
    try:
        with open(input_file, 'r') as csvfile:
            lines = csvfile.readlines()
            if not lines:
                raise ValueError("The input file is empty.")
            # Check if the first line is a header
            timestamps = []
            hd_free_bytes = []
            ram_used_bytes = []
            for line in lines:
                cols = line.strip().split(';')
                if len(cols) != 4:
                    break
                timestamp = cols[0] + ' ' + cols[1]
                timestamps.append(datetime.datetime.strptime(timestamp, '%Y-%m-%d %H:%M:%S'))
                hd_free_bytes.append(int(cols[2]))
                ram_used_bytes.append(int(cols[3]))
        

        hd_free_bytes = np.array(hd_free_bytes)
        ram_used_bytes = np.array(ram_used_bytes)

        # Convert timestamps to a format suitable for plotting
        plt.figure(figsize=(10, 6))
        plt.subplot(2, 1, 1)
        plt.plot(timestamps, hd_free_bytes / 1024**3, label='HD Free Bytes', color='blue')
        ## Y axis show format .2f
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.2f}'))

        plt.xlabel('Time')
        plt.title('HD Free  Over Time')
        plt.ylabel('Free Space [GB]')
        plt.xticks(rotation=45)
        plt.grid()
        plt.legend()
        plt.subplot(2, 1, 2)
        plt.plot(timestamps, ram_used_bytes / 1024**3, label='RAM Used Bytes', color='green')
        plt.xlabel('Time')
        plt.title('RAM Used Over Time')
        plt.ylabel('Used RAM (GB)')
        plt.xticks(rotation=45)
        plt.grid()
        plt.legend()
        plt.tight_layout()
        if output_file:
            plt.savefig(output_file)
        else:
            plt.show()
    except Exception as e:
        print(show_exc(e))

def main():
    parser = argparse.ArgumentParser(description='Plot monitor data from a CSV file.')
    parser.add_argument('input_file', type=str, help='Path to the input CSV file.')
    parser.add_argument('--output', type=str, help='Path to save the output image file.', default=None)
    parser.add_argument('--folder', type=str, help='Folder to save the output image file.', default=None)
    parser.add_argument('--rowheader', action='store_true', help='Indicate if the first row is a header.')

    

    args = parser.parse_args()
    output_path = args.output
    if args.output is not None and args.folder is not None:
        output_path = f"{args.folder}/{args.output}"
    if args.output is None and args.folder is not None:
        basename = os.path.basename(args.input_file)

        output_path = f"{args.folder}/{basename}.png"

    
    plot_monitor_data(args.input_file, output_path, args.rowheader)
if __name__ == '__main__':
    main()