import os
import time
import csv
import argparse
import psutil
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.animation import FuncAnimation
import threading

def is_process_running(pid):
    """Check if a process with given PID is running"""
    return psutil.pid_exists(pid)

def get_free_space_bytes(path):
    """Return available disk space in bytes for a given path"""
    return psutil.disk_usage(path).free

def write_to_csv(csv_file, date, time_str, free_space):
    """Append a new entry to the CSV file"""
    with open(csv_file, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=';')

        writer.writerow([date, time_str, free_space])



def read_csv_data(csv_file):

    """Read CSV data and return timestamps and free space values in GB"""
    timestamps = []
    free_space_gb = []
    if os.path.exists(csv_file):
        with open(csv_file, newline='') as f:
            reader = csv.reader(f, delimiter=';')
            for row in reader:
                datetime_str = f"{row[0]} {row[1]}"
                dt = datetime.strptime(datetime_str, "%Y-%m-%d %H:%M:%S")
                space_gb = int(row[2]) / (1024 ** 3)
                timestamps.append(dt)
                free_space_gb.append(space_gb)
    return timestamps, free_space_gb

def live_plot(csv_file):
    """Display a live plot with mouse hover tooltip"""
    plt.style.use("seaborn")
    fig, ax = plt.subplots()
    annotation = ax.annotate("", xy=(0, 0), xytext=(15, 15), textcoords="offset points",
                             bbox=dict(boxstyle="round", fc="w"),
                             arrowprops=dict(arrowstyle="->"))
    annotation.set_visible(False)
    x_data, y_data = [], []

    def update(frame):
        nonlocal x_data, y_data
        ax.clear()
        x_data, y_data = read_csv_data(csv_file)
        if x_data:
            ax.plot(x_data, y_data, marker='o', linestyle='-')
            ax.set_xlabel("Date and Time")
            ax.set_ylabel("Free Space (GB)")
            ax.set_title("Free Disk Space Over Time")
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            fig.autofmt_xdate()
            annotation.set_visible(False)

    def on_mouse_move(event):
        if not x_data:
            return
        if event.inaxes != ax:
            annotation.set_visible(False)
            fig.canvas.draw_idle()
            return

        # Find closest point
        mouse_x = mdates.num2date(event.xdata)
        distances = [abs((dt - mouse_x).total_seconds()) for dt in x_data]
        index = distances.index(min(distances))
        closest_x = x_data[index]
        closest_y = y_data[index]

        # Update annotation
        annotation.xy = (closest_x, closest_y)
        text = f"{closest_x.strftime('%Y-%m-%d %H:%M:%S')}\n{closest_y:.2f} GB"
        annotation.set_text(text)
        annotation.set_visible(True)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", on_mouse_move)
    ani = FuncAnimation(fig, update, interval=10000)  # update every 10s
    plt.show()


def generate_default_filename():
    """Generate default filename using current datetime"""
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M")
    return f"{timestamp}_fs.log"

def main():
    parser = argparse.ArgumentParser(description="Monitor disk space while a process is running")
    parser.add_argument("pid", type=int, help="PID of the process to monitor")
    parser.add_argument("path", help="Path to check free disk space")
    parser.add_argument("--csv", help="Optional output CSV file name")
    args = parser.parse_args()

    csv_file = args.csv if args.csv else generate_default_filename()

    print(f"Monitoring process PID: {args.pid}")
    print(f"Monitoring path: {args.path}")
    print(f"Logging to file: {csv_file}")
    if not os.path.exists(args.path):
        print(f"Invalid path: {args.path}")
        return

    # Start live plot in a background thread
    threading.Thread(target=live_plot, args=(csv_file,), daemon=True).start()
    while is_process_running(args.pid):
        now = datetime.now()
        date = now.strftime("%Y-%m-%d")
        time_str = now.strftime("%H:%M:%S")
        free_space = get_free_space_bytes(args.path)
        write_to_csv(csv_file, date, time_str, free_space)
        time.sleep(10)
    print("Process ended. Monitor stopped.")

if __name__ == "__main__":

    main()


