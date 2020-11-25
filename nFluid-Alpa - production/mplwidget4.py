# ------------------------------------------------- -----
# -------------------- mplwidget.py --------------------
# -------------------------------------------------- ----
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib.figure import Figure

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as
                                                NavigationToolbar)


class MplWidget4(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.canvas = FigureCanvas(Figure())
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.canvas.Axes3D = self.canvas.figure.add_subplot(221,projection='3d')
        self.canvas.Axes3D2 = self.canvas.figure.add_subplot(222,projection='3d')
        self.canvas.Axes3D3 = self.canvas.figure.add_subplot(223,projection='3d')
        self.canvas.Axes3D4 = self.canvas.figure.add_subplot(224,projection='3d')
        self.canvas.figure.tight_layout(pad=2.0)
        self.canvas.figure.subplots_adjust(left=0, right=0.98, top=1, bottom=0.058, wspace=0.140, hspace=0.140)
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        vertical_layout.addWidget((self.toolbar))
        self.setLayout(vertical_layout)


