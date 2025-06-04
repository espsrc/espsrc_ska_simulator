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
    return psutil.pid_exists(pid)

def get_free_disk_bytes(path):
    return psutil.disk_usage(path).free

def get_used_ram_bytes():
    return psutil.virtual_memory().used

def write_to_csv(csv_file, date, time_str, free_disk, used_ram):
    with open(csv_file, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow([date, time_str, free_disk, used_ram])

def read_csv_data(csv_file):
    timestamps, disk_gb, ram_gb = [], [], []
    if os.path.exists(csv_file):
        with open(csv_file, newline='') as f:
            reader = csv.reader(f, delimiter=';')
            for row in reader:
                dt = datetime.strptime(f"{row[0]} {row[1]}", "%Y-%m-%d %H:%M:%S")
                disk = int(row[2]) / (1024**3)
                ram = int(row[3]) / (1024**3)
                timestamps.append(dt)
                disk_gb.append(disk)
                ram_gb.append(ram)
    return timestamps, disk_gb, ram_gb

def generate_default_filename():
    now = datetime.now()
    return f"{now.strftime('%Y%m%d_%H%M')}_monitor.log"

def live_plot(csv_file):
    #plt.style.use("seaborn")
    fig, (ax_ram, ax_disk) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.tight_layout(pad=3.0)

    annotation_ram = ax_ram.annotate("", xy=(0, 0), xytext=(15, 15), textcoords="offset points",
                                     bbox=dict(boxstyle="round", fc="w"),
                                     arrowprops=dict(arrowstyle="->"))
    annotation_disk = ax_disk.annotate("", xy=(0, 0), xytext=(15, 15), textcoords="offset points",
                                       bbox=dict(boxstyle="round", fc="w"),
                                       arrowprops=dict(arrowstyle="->"))
    annotation_ram.set_visible(False)
    annotation_disk.set_visible(False)

    x_data, y_disk, y_ram = [], [], []

    def update(frame):
        nonlocal x_data, y_disk, y_ram
        ax_ram.clear()
        ax_disk.clear()
        x_data, y_disk, y_ram = read_csv_data(csv_file)

        if x_data:
            ax_ram.plot(x_data, y_ram, marker='o')
            ax_ram.set_ylabel("Used RAM (GB)")
            ax_ram.set_title("Used RAM Over Time")

            ax_disk.plot(x_data, y_disk, marker='o')
            ax_disk.set_ylabel("Free Disk (GB)")
            ax_disk.set_title("Available Disk Space Over Time")

            ax_disk.set_xlabel("Time")
            ax_disk.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            fig.autofmt_xdate()

            annotation_ram.set_visible(False)
            annotation_disk.set_visible(False)

    def on_mouse_move(event):
        if not x_data or event.inaxes not in [ax_ram, ax_disk]:
            annotation_ram.set_visible(False)
            annotation_disk.set_visible(False)
            fig.canvas.draw_idle()
            return

        mouse_time = mdates.num2date(event.xdata)
        diffs = [abs((dt - mouse_time).total_seconds()) for dt in x_data]
        idx = diffs.index(min(diffs))

        dt = x_data[idx]
        ram_val = y_ram[idx]
        disk_val = y_disk[idx]

        if event.inaxes == ax_ram:
            annotation = annotation_ram
            value = ram_val
            label = "RAM"
        else:
            annotation = annotation_disk
            value = disk_val
            label = "Disk"

        annotation.xy = (dt, value)
        annotation.set_text(f"{dt.strftime('%Y-%m-%d %H:%M:%S')}\n{label}: {value:.2f} GB")
        annotation.set_visible(True)
        fig.canvas.draw_idle()

    #fig.canvas.mpl_connect("motion_notify_event", on_mouse_move)
    ani = FuncAnimation(fig, update, interval=10000)
    plt.show()

def monitor_loop(pid, path, csv_file):
    """Background loop that logs RAM and disk usage"""
    while is_process_running(pid):
        now = datetime.now()
        date = now.strftime("%Y-%m-%d")
        time_str = now.strftime("%H:%M:%S")
        free_disk = get_free_disk_bytes(path)
        used_ram = get_used_ram_bytes()
        write_to_csv(csv_file, date, time_str, free_disk, used_ram)
        time.sleep(10)
    print("Process finished. Monitoring stopped.")

def main():
    parser = argparse.ArgumentParser(description="Monitor free RAM and disk space while a process runs")
    parser.add_argument("pid", type=int, help="PID of the process to monitor")
    parser.add_argument("--path", help="Path to check disk space")
    parser.add_argument("--csv", help="CSV output file name")
    args = parser.parse_args()

    csv_file = args.csv if args.csv else generate_default_filename()
    path = args.path if args.path else os.getcwd()

    print(f"Monitoring PID: {args.pid}")
    print(f"Monitoring path: {path}")
    print(f"Logging to: {csv_file}")

    if not os.path.exists(path):
        print(f"Invalid path: {path}")
        return
    
    monitoring_thread = threading.Thread(target=monitor_loop, args=(args.pid, path, csv_file), daemon=True)
    monitoring_thread.start()
    # 👇 Run the plot in the main thread
    live_plot(csv_file)
    # Read CSV data and save the plot into a file
    timestamps, disk_gb, ram_gb = read_csv_data(csv_file)
    fig, (ax_ram, ax_disk) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.tight_layout(pad=3.0)
    ax_ram.plot(timestamps, ram_gb, marker='o')
    ax_ram.set_ylabel("Used RAM (GB)")
    ax_ram.set_title("Used RAM Over Time")
    ax_disk.plot(timestamps, disk_gb, marker='o')
    ax_disk.set_ylabel("Free Disk (GB)")
    ax_disk.set_title("Available Disk Space Over Time")
    ax_disk.set_xlabel("Time")
    ax_disk.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    fig.autofmt_xdate()
    plt.savefig(f"{csv_file}.png")
    print(f"Plot saved as {csv_file}.png")

if __name__ == "__main__":
    main()
