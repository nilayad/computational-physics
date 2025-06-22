import matplotlib.pyplot as plt

def tenplot(m_state):
    
# Flatten all (tensor index, slice index) pairs into one list
    all_slices = []
    for j, tensor in m_state.items():
        num_slices = tensor.shape[1]  # Only 2 slices per tensor
        for i in range(num_slices):
            all_slices.append((j, i, tensor[:, i, :]))

# Plot in rows of 6
    num_images = len(all_slices)
    cols = 6
    rows = (num_images + cols - 1) // cols  # ceiling division

    fig, axs = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))

# Ensure axs is always 2D
    axs = axs.reshape((rows, cols)) if rows > 1 else [axs]

    for idx, (j, i, img) in enumerate(all_slices):
        row, col = divmod(idx, cols)
        ax = axs[row][col]
        im = ax.imshow(img, aspect='auto')
        ax.set_title(f"Tensor {j}, Slice {i}")
    
    #fig.colorbar(im)
# Hide any unused subplots
    for idx in range(len(all_slices), rows * cols):
        row, col = divmod(idx, cols)
        axs[row][col].axis('off')
    cbar = fig.colorbar(im, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('Intensity')
    plt.grid()
    plt.tight_layout()
    plt.show()


import matplotlib.pyplot as plt
import numpy as np

def tenplot2(m_state):
    # Flatten all (tensor index, slice index) pairs into one list
    all_slices = []
    for j, tensor in m_state.items():
        num_slices =  tensor.shape[1]  # Only 2 slices per tensor
        for i in range(num_slices):
            all_slices.append((j, i, tensor[:, i, :]))

    # Plot in rows of 6
    num_images = len(all_slices)
    cols = 6
    rows = (num_images + cols - 1) // cols  # ceiling division

    fig, axs = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))

    # Ensure axs is 2D
    if rows == 1:
        axs = axs[np.newaxis, :]

    # Find global min and max across all images for consistent color scaling
    #all_data = np.array([img for (_, _, img) in all_slices])
    vmin = min(img.min() for _, _, img in all_slices)
    vmax = max(img.max() for _, _, img in all_slices)


    im = None
    for idx, (j, i, img) in enumerate(all_slices):
        row, col = divmod(idx, cols)
        ax = axs[row][col]
        im = ax.imshow(img, aspect='auto', vmin=vmin, vmax=vmax, cmap='viridis')
        ax.set_title(f"Tensor {j}, Slice {i}")

        h, w = img.shape
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width_px = bbox.width * fig.dpi
        height_px = bbox.height * fig.dpi

        cell_width = width_px / w
        cell_height = height_px / h

        # Choose font size to be ~80% of the smaller dimension of a cell
        cell_fontsize = 0.2 * min(cell_width, cell_height)
        
        # Add text annotations (rounded to 2 decimals)
        for (m, n), value in np.ndenumerate(img):
            ax.text(n, m, f"{value:.1f}", ha='center', va='center',
                    fontsize=cell_fontsize, color='white' if value < (vmin + vmax) / 2 else 'black')

    # Hide any unused subplots
    for idx in range(len(all_slices), rows * cols):
        row, col = divmod(idx, cols)
        axs[row][col].axis('off')

    # Create a single colorbar for the entire figure
    cbar = fig.colorbar(im, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('Value')

    #plt.tight_layout()
    plt.show()
