import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class DetailedPlots():

    def __init__(self, model):
        self.model = model

    def plot_aerogrid(self, scalars=None, colormap='plasma'):
        """
        Plot the aerodynamic grid using Matplotlib's 3D plotting of polygons.

        Parameters:
        -----------
        scalars : array-like, optional
            Values to use for coloring the panels
        colormap : str, default='plasma'
            Matplotlib colormap name
        """
        # Create figure and 3D axes
        fig = plt.figure(figsize=(16, 9))
        ax = fig.add_subplot(111, projection='3d')

        # Get points and create panels
        points = self.model.aerogrid['cornerpoint_grids'][:, (1, 2, 3)]
        panels = []
        for panel in self.model.aerogrid['cornerpoint_panels']:
            # Get corner points for each panel
            panel_points = [points[np.where(self.model.aerogrid['cornerpoint_grids'][:, 0] == id)[0][0]]
                            for id in panel]
            panels.append(panel_points)

        # Create 3D collection of polygons
        poly3d = Poly3DCollection(panels, edgecolor='black', linewidth=0.2)

        # Set face colors based on scalars if provided
        if scalars is not None:
            cmap = plt.get_cmap(colormap)
            norm = plt.Normalize(vmin=np.min(scalars), vmax=np.max(scalars))
            poly3d.set_facecolor(cmap(norm(scalars)))
            poly3d.set_alpha(1.0)

            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax)
        else:
            poly3d.set_facecolor('white')
            poly3d.set_alpha(0.7)

        # Add collection to axes
        ax.add_collection3d(poly3d)

        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Set equal aspect ratio for all axes
        # This needs to be calculated manually with the center as reference point
        x_range = points[:, 0].max() - points[:, 0].min()
        y_range = points[:, 1].max() - points[:, 1].min()
        z_range = points[:, 2].max() - points[:, 2].min()
        max_range = max(x_range, y_range, z_range)

        mid_x = (points[:, 0].max() + points[:, 0].min()) * 0.5
        mid_y = (points[:, 1].max() + points[:, 1].min()) * 0.5
        mid_z = (points[:, 2].max() + points[:, 2].min()) * 0.5

        ax.set_xlim(mid_x - max_range * 0.5, mid_x + max_range * 0.5)
        ax.set_ylim(mid_y - max_range * 0.5, mid_y + max_range * 0.5)
        ax.set_zlim(mid_z - max_range * 0.5, mid_z + max_range * 0.5)

        # Set equal box aspect
        ax.set_box_aspect((1, 1, 1))

        # Set default view angle
        ax.view_init(elev=40, azim=-120)

        plt.show()
