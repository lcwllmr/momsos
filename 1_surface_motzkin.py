import numpy as np
import matplotlib.pyplot as plt
from helper import motzkin_xy

def main():
    x = np.linspace(-1.5, 1.5, 401)
    y = np.linspace(-1.5, 1.5, 401)
    X, Y = np.meshgrid(x, y)
    Z = motzkin_xy(X, Y)

    # Focus on the interesting range around the near-zero valleys.
    Zc = np.clip(Z, -0.2, 3.0)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Zc, rstride=4, cstride=4, linewidth=0, antialiased=True)
    ax.set_title("Motzkin polynomial (clipped for visibility)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("m(x,y)")
    ax.view_init(elev=28, azim=-55)

    plt.tight_layout()
    plt.savefig("motzkin_surface.png", dpi=200)
    print("Wrote motzkin_surface.png")

if __name__ == "__main__":
    main()

