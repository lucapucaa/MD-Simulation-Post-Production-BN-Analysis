import argparse
import sys
import matplotlib
# Force headless image generation before importing pyplot
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
import numpy as np


def plot_xvg():
    # 1. Set up command-line arguments
    parser = argparse.ArgumentParser(
        description="Plot any GROMACS .xvg file dynamically, supporting multi-column data."
    )
    parser.add_argument(
        "-f", "--file", required=True, help="Path to the input .xvg file"
    )
    parser.add_argument(
        "-t", "--title", default="GROMACS Analysis", help="Title of the plot"
    )
    parser.add_argument("-x", "--xlabel", default="Time (ns)", help="X-axis label")
    parser.add_argument("-y", "--ylabel", default="Value (nm)", help="Y-axis label")
    parser.add_argument(
        "-l", "--label", default="Data", 
        help="Comma-separated legend labels for data lines (e.g., 'Backbone,All-Atom')"
    )
    parser.add_argument(
        "-o", "--output", default="plot.png", help="Output image file name"
    )
    parser.add_argument(
        "--ps", action="store_true", 
        help="Keep X-axis in picoseconds (default converts raw GROMACS ps to ns)"
    )

    args = parser.parse_args()

    # 2. Load the data safely, bypassing GROMACS header annotations
    try:
        data = np.loadtxt(args.file, comments=["@", "#"])
    except Exception as e:
        print(f"Error reading file '{args.file}': {e}")
        sys.exit(1)

    # Validate that data is not completely empty
    if data.size == 0 or data.ndim == 0:
        print(f"Error: File '{args.file}' contains no parseable numeric matrix data.")
        sys.exit(1)

    # 3. Handle X-axis and unit transformation
    # If the array is flat (1D), treat it as a single data column map
    x = data if data.ndim == 1 else data[:, 0]
    
    # Auto-convert time arrays from ps to ns unless overridden by the user
    if not args.ps and "time" in args.xlabel.lower() and data.ndim > 1:
        x = x / 1000.0  

    plt.figure(figsize=(8, 5))

    # 4. Handle Multi-Column Data Dynamically
    labels = [lbl.strip() for lbl in args.label.split(",")]

    if data.ndim > 1 and data.shape[1] > 2:
        # Multi-column array routing (e.g., Radius of Gyration decomposition axes)
        for i in range(1, data.shape[1]):
            lbl = labels[i-1] if i-1 < len(labels) else f"Column {i}"
            plt.plot(x, data[:, i], label=lbl, linewidth=1.5)
    elif data.ndim > 1 and data.shape[1] == 2:
        # Standard single Y-column data (e.g., simple Backbone RMSD or RMSF profiles)
        plt.plot(x, data[:, 1], label=labels[0], color="blue", linewidth=1.5)
    else:
        # Fallback track for simple array lists
        plt.plot(x, x, label=labels[0], color="blue", linewidth=1.5)

    # Apply variable titles and labels passed from the command line interface
    plt.title(args.title, fontsize=13, fontweight='bold')
    plt.xlabel(args.xlabel, fontsize=11)
    plt.ylabel(args.ylabel, fontsize=11)

    # 5. Styling, cleanup, and file output
    plt.legend(loc="best", frameon=True, facecolor="white", edgecolor="none")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    plt.savefig(args.output, dpi=300)
    print(f"Success! Plot saved as '{args.output}'")
    plt.close()  # Clear memory allocations for automated pipeline loops


if __name__ == "__main__":
    plot_xvg()
