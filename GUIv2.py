from __future__ import print_function, division, absolute_import
from frankenScope_well_reshape import filenamesToDict
import sys
from PyQt4 import QtGui, QtCore


class AppWindow(QtGui.QWidget):

    def __init__(self):
        super(AppWindow, self).__init__()
        self.indir = None
        self.outdir = None
        self.initUI()

    def initUI(self):
        #All visible lables in order top->bottom
        lbl1 = QtGui.QLabel("Load files from this directory:", self)
        btn1_label = "Select input directory"
        self.indir_lbl = QtGui.QLabel(str(self.indir), self)

        lbl2 = QtGui.QLabel("Save stitched files to this directory:", self)
        btn2_label = "Select output directory"
        self.outdir_lbl = QtGui.QLabel(str(self.outdir), self)

        #Tooltips for the buttons
        btn1_tooltip = "Select the directory where yor OME-TIFF files are " \
                       "stored"
        btn2_tooltip = "Select the directory where your stitched files will " \
                       "be stored"

        #Make buttons work
        btn1 = QtGui.QPushButton(btn1_label, self)
        btn1.clicked.connect(self.select_indir)
        btn1.setToolTip(btn1_tooltip)
        btn1.resize(btn1.sizeHint())


        btn2 = QtGui.QPushButton(btn2_label, self)
        btn2.clicked.connect(self.select_outdir)
        btn2.setToolTip(btn1_tooltip)
        btn2.resize(btn2.sizeHint())


        #QtGui.QToolTip.setFont(QtGui.QFont("SansSerif", 12))
        #self.setToolTip("This is a work in <i>progress</i>...")

        #exitAction = QtGui.QAction(QtGui.QIcon("logo.png"), "&Stop", self)
        #exitAction.setShortcut("Ctrl+Q")
        #exitAction.setStatusTip("Exit application")
        #exitAction.triggered.connect(QtGui.qApp.quit)

        #self.statusBar().showMessage("Ready")

        #menubar = self.menuBar()
        #fileMenu = menubar.addMenu("&File")
        #fileMenu.addAction(exitAction)



        qbtn = QtGui.QPushButton("Quit", self)
        qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        qbtn.resize(qbtn.sizeHint())
        #qbtn.move(100, 50)

        runButton = QtGui.QPushButton("Run...", self)
        runButton.clicked.connect(self.center)
        qbtn.resize(qbtn.sizeHint())

        #self.dialogbtn = QtGui.QPushButton("Dialog", self)
        #self.dialogbtn.clicked.connect(self.showDialog)
        #self.le = QtGui.QLineEdit(self)

        #grid = QtGui.QGridLayout()
        #self.setLayout(grid)
        #names =["Ass", "Hat", '3', '4', '',
        #        "Kikkeli", "bajs", "snopp", "666", "667"]

        #positions = [(i, j) for i in range(2) for j in range(5)]

        #for position, name in zip(positions, names):
         #   if name == '':
          #      continue
           # button = QtGui.QPushButton(name)
            #grid.addWidget(button, *position)
        runQuitBox = QtGui.QHBoxLayout()
        runQuitBox.addWidget(runButton)
        runQuitBox.addWidget(qbtn)

        verticalWidgets = [lbl1, self.indir_lbl, btn1, lbl2, self.outdir_lbl,
                           btn2]
        vbox = QtGui.QVBoxLayout()
        for widget in verticalWidgets:
            vbox.addWidget(widget)

        vbox.addStretch(1)
        vbox.addLayout(runQuitBox)
        #vbox = QtGui.QVBoxLayout()
        #vbox.addStretch(1)


        self.setLayout(vbox)
        self.setGeometry(300, 300, 300, 500)
        self.setWindowTitle("Stitch buddy")
        self.setWindowIcon(QtGui.QIcon("logo.png"))
        #self.center()
        self.show()




        self.show()

    def showDialog(self):

        col = QtGui.QColorDialog.getColor()

        text, ok = QtGui.QInputDialog.getText(self, "input dialog",
                                              "Enter name:")
        if ok:
            self.le.setText(str(text))


    def select_indir(self):
        self.indir = QtGui.QFileDialog.getExistingDirectory(None,
                                                            'Select input folder...','',
                                                            QtGui.QFileDialog.ShowDirsOnly)

        print(self.indir)
        self.indir_lbl.setText(str(self.indir))

    def select_outdir(self):
        self.outdir = QtGui.QFileDialog.getExistingDirectory(None,
                                                            'Select output folder...',
                                                            '',
                                                            QtGui.QFileDialog.ShowDirsOnly)

        self.outdir_lbl.setText(str(self.outdir))

    def getWelldict(self):

        self.wellDict = filenamesToDict(self.indir)

    def center(self):

        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def closeEvent(self, event):

        reply = QtGui.QMessageBox.question(self, "Message",
                                           "Are you sure you want to quit?",
                                           QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                                           QtGui.QMessageBox.Yes)
        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


def main():

    app = QtGui.QApplication(sys.argv)
    showMe = AppWindow()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
