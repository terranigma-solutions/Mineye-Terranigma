import numpy as np
from bayseg import BaySeg
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import time

# Load your prepared stack
print("Loading data...")
layers = np.load("sentinel2_bayseg_input.npy")
print(f"Input data shape: {layers.shape}")

# Initialize segmenter
n_classes = 6
print(f"Running segmentation with {n_classes} classes...")
start_time = time.time()
seg = BaySeg(data=layers, n_labels=n_classes, beta_init=10)

# Fit the model
labels = seg.fit(n=200, beta_jump_length=0.1)
print(f"Segmentation completed in {time.time() - start_time:.2f} seconds")

# Save result
np.save("bayseg_lithology_labels.npy", labels)
print("Results saved to bayseg_lithology_labels.npy")

# Plot diagnostic information
seg.diagnostics_plot(transpose=False)

# Plot segmentation results
plt.figure(figsize=(10, 8))
plt.title(f"Segmentation Result ({n_classes} classes)")
cmap = ListedColormap(plt.cm.tab10(range(n_classes)))
plt.imshow(labels, cmap=cmap)
plt.colorbar(ticks=range(n_classes), label="Class")
plt.axis('off')

# Plot acceptance ratios
seg.plot_acc_ratios()

# This is required to display plots in PyCharm
plt.tight_layout()
plt.show()