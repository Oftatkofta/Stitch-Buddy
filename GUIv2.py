from __future__ import print_function, division, absolute_import
from frankenScope_well_reshape import filenamesToDict
import sys
from PyQt4 import QtGui, QtCore
import re

class EmittingStream(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, text):
        self.textWritten.emit(str(text))

class ErrorEmittingStream(QtCore.QObject):

    textWritten = QtCore.pyqtSignal(str)

    def write(self, errortext):
        self.textWritten.emit(str(errortext))



class AppWindow(QtGui.QWidget):

    def __init__(self):
        super(AppWindow, self).__init__()

        self.indir = None
        self.outdir = None
        self.wellDict = None

        # Installs custom output streams
        sys.stdout = EmittingStream(textWritten=self.normalOutputWritten)
        sys.stderr = ErrorEmittingStream(textWritten=self.errorOutputWritten)

        self.initUI()

    def __del__(self):
        # Restore sys.stdout and sys.stderr
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__


    def initUI(self):

        #All visible lables in order top->bottom
        lbl1 = QtGui.QLabel("Load files from this directory:", self)
        btn1_label = "Select input directory"
        self.indir_lbl = QtGui.QLabel(str(self.indir), self)

        lbl2 = QtGui.QLabel("Save stitched files to this directory:", self)
        btn2_label = "Select output directory"
        self.outdir_lbl = QtGui.QLabel(str(self.outdir), self)

        btn3_label = "Rename wells (optional)"

        #Tooltips for the buttons
        btn1_tooltip = "Select the directory where yor OME-TIFF files are " \
                       "stored"
        btn2_tooltip = "Select the directory where your stitched files will " \
                       "be stored"
        btn3_tooltip = "Opens an editor where you can write the names of the" \
                       " stitched wells"

        #Make buttons work
        btn1 = QtGui.QPushButton(btn1_label, self)
        btn1.clicked.connect(self.select_indir)
        btn1.setToolTip(btn1_tooltip)
        btn1.resize(btn1.sizeHint())

        btn2 = QtGui.QPushButton(btn2_label, self)
        btn2.clicked.connect(self.select_outdir)
        btn2.setToolTip(btn1_tooltip)
        btn2.resize(btn2.sizeHint())

        btn3 = QtGui.QPushButton(btn3_label, self)
        btn3.clicked.connect(self.rename_wells)
        btn3.setToolTip(btn3_tooltip)
        btn3.resize(btn3.sizeHint())

        #textbox to log output and errors from component functions
        self.logOutput = QtGui.QTextEdit()
        self.logOutput.setReadOnly(True)


        qbtn = QtGui.QPushButton("Quit", self)
        qbtn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        qbtn.resize(qbtn.sizeHint())

        runButton = QtGui.QPushButton("Run...", self)
        runButton.clicked.connect(self.center)
        qbtn.resize(qbtn.sizeHint())

        #Group run and quit buttons into one row
        runQuitBox = QtGui.QHBoxLayout()
        runQuitBox.addWidget(runButton)
        runQuitBox.addWidget(qbtn)

        #Separator line
        line1 = QtGui.QFrame()
        line1.setFrameShape(QtGui.QFrame.HLine)
        line1.setFrameShadow(QtGui.QFrame.Sunken)


        verticalWidgets = [lbl1, self.indir_lbl, btn1, lbl2, self.outdir_lbl,
                           btn2, line1, btn3]

        vbox = QtGui.QVBoxLayout()
        for widget in verticalWidgets:
            vbox.addWidget(widget)

        vbox.addWidget(self.logOutput)
        vbox.addStretch(1)
        vbox.addLayout(runQuitBox)


        self.setLayout(vbox)
        self.setGeometry(300, 300, 300, 500)
        self.setWindowTitle("Stitch buddy")
        self.setWindowIcon(QtGui.QIcon("logo.png"))
        #self.center()
        self.show()

    def normalOutputWritten(self, text):
        """Append text to the QTextEdit."""

        self.logOutput.insertPlainText(text)
        sb= self.logOutput.verticalScrollBar()
        sb.setValue(sb.maximum())
        #self.logOutput.setTextCursor(cursor)
        #self.logOutput.ensureCursorVisible()

    def errorOutputWritten(self, errortext):
        """Append red error text to the QTextEdit."""

        #sets fontcolor to red for warnings
        color = QtGui.QColor(255, 0, 0)
        self.logOutput.setTextColor(color)

        #Write ouput to log
        self.logOutput.insertPlainText(errortext)

        #Set fontcolor back to black
        color = QtGui.QColor(0, 0, 0)
        self.logOutput.setTextColor(color)

        #Autoscroll the text
        sb= self.logOutput.verticalScrollBar()
        sb.setValue(sb.maximum())


    def showDialog(self):

        col = QtGui.QColorDialog.getColor()

        text, ok = QtGui.QInputDialog.getText(self, "input dialog",
                                              "Enter name:")
        if ok:
            self.le.setText(str(text))



    def select_indir(self):
        self.indir = QtGui.QFileDialog.getExistingDirectory(None,
                                                            'Select input folder...',
                                                            '/Volumes/HDD/Huygens_SYNC/Raw OME files for test/2Z_2C_2x2gridx2_2T_2',
                                                            QtGui.QFileDialog.ShowDirsOnly)

        self.wellDict = filenamesToDict(str(self.indir))
        self.indir_lbl.setText(str(self.indir))
        boldFont = QtGui.QFont()
        boldFont.setBold(True)
        self.indir_lbl.setFont(boldFont)

    def select_outdir(self):
        self.outdir = QtGui.QFileDialog.getExistingDirectory(None,
                                                            'Select output folder...',
                                                            '',
                                                            QtGui.QFileDialog.ShowDirsOnly)

        self.outdir_lbl.setText(str(self.outdir))
        boldFont = QtGui.QFont()
        boldFont.setBold(True)
        self.outdir_lbl.setFont(boldFont)


    def rename_wells(self):

        if self.wellDict == None:
            print("No indir selected!")

        else:
            for key in self.wellDict:
                print(self.wellDict[key])


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
