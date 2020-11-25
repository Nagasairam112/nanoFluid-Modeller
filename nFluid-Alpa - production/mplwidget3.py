# ------------------------------------------------- -----
# -------------------- mplwidget.py --------------------
# -------------------------------------------------- ----
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as
                                                NavigationToolbar)


class MplWidget3(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.canvas = FigureCanvas(Figure())
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.canvas.axes = self.canvas.figure.add_subplot(221)
        self.canvas.axes2 = self.canvas.figure.add_subplot(222)
        self.canvas.axes3 = self.canvas.figure.add_subplot(223)
        self.canvas.axes4 = self.canvas.figure.add_subplot(224)
        self.canvas.figure.tight_layout(pad=2.0)

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        vertical_layout.addWidget((self.toolbar))
        self.setLayout(vertical_layout)


