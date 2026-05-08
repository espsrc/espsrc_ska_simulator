import csv
import os
from datetime import datetime
from io import StringIO
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import argparse

def read_csv(csv_file):
    timestamps = []
    memory = []
    if not os.path.exists(csv_file):
        return timestamps, memory
    with open(csv_file, newline='') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            if len(row) != 3:
                continue
            try:
                dt = datetime.strptime(f"{row[0].strip()} {row[1].strip()}", "%Y-%m-%d %H:%M:%S")
                mem = float(row[2].strip())
                timestamps.append(dt)
                memory.append(mem)
            except Exception:
                continue
    return np.array(timestamps), np.array(memory)



def update_plot(csv_file, params):
    [fig, ax, line] = params
    x, y = read_csv(csv_file)
    if x.all() and y.all():
        y /= (1024.**2)
        line.set_data(x, y)
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

    return [fig, ax, line]

# Mouse move event to show coordinates
def on_mouse_move2(event):
    if event.xdata and event.ydata:
        dt_str = mdates.num2date(event.xdata).strftime('%Y-%m-%d %H:%M:%S')
        mem_str = f"{event.ydata:.2f} GB"
        coord_text.set_text(f"Datetime: {dt_str}\nMemory: {mem_str}")
        fig.canvas.draw_idle()

# Mouse move event to show coordinates
def on_mouse_move(event):
    if event.xdata and event.ydata:
        ax = event.inaxes
        dt_str = mdates.num2date(event.xdata).strftime('%Y-%m-%d %H:%M:%S')
        mem_str = f"{event.ydata:.2f} GB"
        coord_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=10,
                 verticalalignment='top', bbox=dict(boxstyle="round", facecolor='wheat', alpha=0.5))
        coord_text.set_text(f"Datetime: {dt_str}\nMemory: {mem_str}")
        event.canvas.draw_idle()


def main():
    parser = argparse.ArgumentParser(description="Monitor and plot memory usage from a CSV file.")
    parser.add_argument("csv_file", type=str, help="Path to the CSV file containing memory usage data.")
    parser.add_argument("--sep", default=";", help="CSV separator (default: ';')")
    parser.add_argument("--interval", type=float, default=30.0, help="Refresh interval in seconds.")
    args = parser.parse_args()

    csv_file = args.csv_file
# Check if the CSV file exists
    if not os.path.exists(csv_file):
        print(f"Error: The file {csv_file} does not exist.")
        exit(1)

# Setup the plot
    plt.ion()
    fig, ax = plt.subplots()
    line, = ax.plot([], [], '-.')
    ax.set_title("Memory Usage Over Time")
    ax.set_xlabel("Date and Time")
    ax.set_ylabel("Memory (GB)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    fig.autofmt_xdate()
    params = [fig, ax, line]



# Text box to show coordinates

    fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)

    last_mtime = 0

    print ("Click Ctrl + C to finish....")
    while True:
        try:
            mtime = os.path.getmtime(csv_file)
            if mtime != last_mtime:
                last_mtime = mtime
                params = update_plot(csv_file, params)
            plt.pause(args.interval)
        except FileNotFoundError:
            print("CSV file not found. Waiting...")
            time.sleep(args.interval)
        except KeyboardInterrupt:
            print("Exiting...")
            break


if __name__ == "__main__":
    main()
