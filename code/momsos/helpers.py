import sys
import time
from matplotlib import pyplot as plt

def build_plot_in_context(name, size, plot_fn):
    """
    This is just a wrapper function for generating plots that will both
    be viewed immediately for developing and saved in dark and light mode
    theme for publishing with my custom markdown processors.
    Probably useless for other people.
    """
    is_scimd = '--scimd' in sys.argv
    figure = plt.figure(figsize=size)

    if is_scimd:
        plt.style.use("default")
        plot_fn(figure)
        plt.tight_layout()
        plt.savefig(f"{name}.light.png", transparent=True)
        print(f"Wrote {name}.light.png")

        plt.style.use("dark_background")
        plot_fn(figure)
        plt.tight_layout()
        plt.savefig(f"{name}.dark.png", transparent=True)
        print(f"Wrote {name}.dark.png")
    else:
        plt.style.use("default")
        plot_fn(figure)
        plt.tight_layout()
        plt.show()
