import numpy as np
from matplotlib import pyplot as plt
from momsos.helpers import build_plot_in_context
from momsos.poly import Polynomial

def main():
    res = 401
    x = np.linspace(-1.4, 1.4, res)
    y = np.linspace(-1.4, 1.4, res)
    X, Y = np.meshgrid(x, y)

    Z = Polynomial.motzkin()(np.vstack([X.flatten(), Y.flatten()]))
    Z = Z.reshape((res, res))
    Z[(Z < 0) | (Z > 1.2)] = np.nan

    def plot(figure: plt.Figure):
        ax = figure.add_subplot(111, projection="3d")
        ax.plot_surface(X, Y, Z, rcount=100, ccount=100, cmap='jet', antialiased=True)
        ax.set_zlim(0.0, 1.2)
        ax.view_init(elev=28, azim=-55)
    build_plot_in_context("motzkin-surface", (9, 7), plot)


if __name__ == "__main__":
    main()

