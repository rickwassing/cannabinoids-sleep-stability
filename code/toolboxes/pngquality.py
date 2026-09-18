import os
import csv
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# --- Configuration ---
IMG_FOLDER = "/Volumes/sleep/Sleep/3. ACTIVE STUDIES/CUPID/Arousal paper backup_CANSLEEP/inspect/n2aros"   # change this
CSV_FILE = "/Volumes/sleep/Sleep/3. ACTIVE STUDIES/CUPID/Arousal paper backup_CANSLEEP/inspect/selectedaros.csv"

# --- Load previous labels (if any) ---
labels = {}
if os.path.exists(CSV_FILE):
    with open(CSV_FILE, newline="") as f:
        reader = csv.reader(f)
        labels = {row[0]: row[1] for row in reader}

# --- Get list of images ---
files = sorted([f for f in os.listdir(IMG_FOLDER) if f.lower().endswith(".png")])
files_to_do = [f for f in files if f not in labels]

print(f"{len(files)} images total, {len(files_to_do)} left to label.")

# --- Set up matplotlib window ---
fig, ax = plt.subplots()
plt.axis("off")

def on_key(event):
    global idx
    if event.key in ["g", "b"]:  # good / bad
        labels[files_to_do[idx]] = event.key
        with open(CSV_FILE, "w", newline="") as f:
            writer = csv.writer(f)
            for k, v in labels.items():
                writer.writerow([k, v])
        print(f"{files_to_do[idx]} → {event.key.upper()}")
        idx += 1
        show_next()
    elif event.key == "q":  # quit
        print("Quitting.")
        plt.close()

def show_next():
    if idx < len(files_to_do):
        ax.clear()
        ax.axis("off")
        img = mpimg.imread(os.path.join(IMG_FOLDER, files_to_do[idx]))
        ax.imshow(img)
        ax.set_title(f"{files_to_do[idx]}  ({idx+1}/{len(files_to_do)})\n[g]=good, [b]=bad, [q]=quit")
        fig.canvas.draw()
    else:
        print("All done!")
        plt.close()

# --- Run ---
idx = 0
fig.canvas.mpl_connect("key_press_event", on_key)
show_next()
plt.show()