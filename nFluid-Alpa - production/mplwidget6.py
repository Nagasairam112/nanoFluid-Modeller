# ------------------------------------------------- -----
# -------------------- mplwidget.py --------------------
# -------------------------------------------------- ----
import matplotlib
from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure


from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as
                                                NavigationToolbar)
from qtconsole.qt import QtCore


class MplWidget6(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.canvas = FigureCanvas(Figure())
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.canvas.A = self.canvas.figure.add_subplot(221)
        self.canvas.A2 = self.canvas.figure.add_subplot(222)
        self.canvas.A3 = self.canvas.figure.add_subplot(223)
        self.canvas.A4 = self.canvas.figure.add_subplot(224)
        self.canvas.figure.tight_layout(pad=2)
        self.canvas.figure.subplots_adjust(left=0.107, right=0.959, top=0.94, bottom=0.132, wspace=0.39,hspace=0.422)
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        vertical_layout.addWidget((self.toolbar))
        self.toolbar.setOrientation(QtCore.Qt.Horizontal)

        self.setLayout(vertical_layout)


