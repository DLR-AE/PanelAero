import numpy as np

try:
    from mayavi import mlab
    from tvtk.api import tvtk
except Exception:
    pass


class DetailedPlots():

    def __init__(self, model):
        self.model = model

    def plot_aerogrid(self, scalars=None, colormap='plasma', embed_in_notebook=False):
        # create the unstructured grid
        points = self.model.aerogrid['cornerpoint_grids'][:, (1, 2, 3)]
        ug = tvtk.UnstructuredGrid(points=points)
        shells = []
        for shell in self.model.aerogrid['cornerpoint_panels']:
            shells.append([np.where(self.model.aerogrid['cornerpoint_grids'][:, 0] == id)[0][0] for id in shell])
        shell_type = tvtk.Polygon().cell_type
        ug.set_cells(shell_type, shells)
        ug.cell_data.scalars = scalars

        # hand over unstructured grid to mayavi
        if embed_in_notebook:
            mlab.init_notebook('png')
        mlab.figure(bgcolor=(1, 1, 1), size=(600, 450))
        src_aerogrid = mlab.pipeline.add_dataset(ug)

        # determine if suitable scalar data is given
        if scalars is not None:
            surface = mlab.pipeline.surface(src_aerogrid, opacity=1.0, line_width=0.5,
                                            colormap=colormap, vmin=scalars.min(), vmax=scalars.max())
            surface.actor.mapper.scalar_visibility = True

            surface.module_manager.scalar_lut_manager.show_legend = True
            surface.module_manager.scalar_lut_manager.data_name = ''
            surface.module_manager.scalar_lut_manager.label_text_property.color = (0, 0, 0)
            surface.module_manager.scalar_lut_manager.label_text_property.font_family = 'arial'
            surface.module_manager.scalar_lut_manager.label_text_property.bold = False
            surface.module_manager.scalar_lut_manager.label_text_property.italic = False
            surface.module_manager.scalar_lut_manager.number_of_labels = 5
        else:
            surface = mlab.pipeline.surface(src_aerogrid, opacity=1.0, line_width=0.5)
            surface.actor.mapper.scalar_visibility = False
        surface.actor.property.edge_visibility = True
        mlab.view(azimuth=60.0, elevation=-65.0, roll=55.0)

        if embed_in_notebook:
            return surface
        else:
            mlab.show()
