"""Module Containing objects for running molecular visualisation in a Jupyter Notebook cell."""

import os
import tempfile

from ipywidgets import interact
# catch for failed import of mdtraj if C++ compiler missing (in Windows)
try:
    from mdtraj import load
    from mdtraj.html import TrajectoryView, enable_notebook
except ImportError:
    print("Import successful except failed to import mdtraj. The notebook visualiser (AMPALViewer) will not function.",
          "\nEnsure the mdtraj module is installed (http://mdtraj.org/latest/). If you are on Windows, also ensure you\
           have a C++ compiler installed (see https://matthew-brett.github.io/pydagogue/python_msvc.html).")


class AMPALViewer(object):
    def __init__(self, ampal_object):
        """Used to visualise an AMPAL Object in a Jupyter Notebook."""
        enable_notebook()
        self.ampal_object = ampal_object
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as t_outf:
                t_outf.write(ampal_object.pdb.encode())
                t_outf.seek(0)
                self.traj = load(t_outf.name)
        finally:
            os.remove(t_outf.name)
        self.ampal_view = TrajectoryView(self.traj)
        self.ampal_view.secondaryStructure = 'ribbon'
        self.ampal_view.background = 'grey'

    def view(self):
        """Starts molecular viewer"""
        return self.ampal_view

    def control_panel(self):
        """Launches control panel for molecular viewer"""
        primary_structure = ('nothing', 'lines', 'stick', 'ball & stick', 'sphere')
        secondary_structure = ('ribbon', 'strand', 'cylinder & plate', 'C alpha trace', 'nothing')
        # Surface is slow to render
        # surface_representation = ('nothing', 'Van der Waals surface', 'solvent excluded surface',
        #                          'solvent accessible surface', 'molecular surface')
        colour_by = ('spectrum', 'chain', 'secondary structure', 'residue', 'polarity', 'atom')
        background_colour = ('grey', 'black', 'white')
        camera_type = ('perspective', 'orthographic')

        panel = interact(self.set_renderer,
                         primary_structure=primary_structure,
                         secondary_structure=secondary_structure,
                         # surface_representation=surface_representation,
                         colour_by=colour_by,
                         background_colour=background_colour,
                         camera_type=camera_type
                         )
        return panel

    def set_renderer(self, primary_structure, secondary_structure, colour_by, background_colour, camera_type):
        """Sets protein representation values in viewer."""
        self.ampal_view.primaryStructure = primary_structure
        self.ampal_view.secondaryStructure = secondary_structure
        self.ampal_view.colorBy = colour_by
        self.ampal_view.background = background_colour
        self.ampal_view.camera = camera_type
        return


__author__ = "Christopher W. Wood"
