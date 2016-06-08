import sys
from PyQt4 import QtGui, QtCore

_version = 0.1

class Window(QtGui.QMainWindow):



    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("MultiWell Mosaic Stitcher")
        self.setWindowIcon(QtGui.QIcon('logo.png'))

        exitAction = QtGui.QAction("&Swixy the app!", self)
        exitAction.setShortcut("Ctrl+E")
        exitAction.setStatusTip("Quit the app")
        exitAction.triggered.connect(self.close_application)

        extractAction = QtGui.QAction("&A filemenu item...", self)
        extractAction.setShortcut("Ctrl+Q")
        extractAction.setStatusTip('Leave The App')
        extractAction.triggered.connect(self.close_application)


        self.statusBar()

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu("&Grilly")
        fileMenu.addAction(exitAction)
        fileMenu.addAction(extractAction)

        self.completed = 0
        self.home()

    def home(self):
        quit_btn = QtGui.QPushButton("Quit", self)
        quit_btn.clicked.connect(self.close_application)
        quit_btn.resize(quit_btn.sizeHint())
        quit_btn.move(0,75)

        extractAction = QtGui.QAction(QtGui.QIcon('logo.png'), "hovertext", self)
        extractAction.triggered.connect(self.close_application)

        self.toolBar = self.addToolBar("Extraction")
        self.toolBar.addAction(extractAction)
        checkbox = QtGui.QCheckBox("Enlage window", self)
        checkbox.move(0,50)
        checkbox.stateChanged.connect(self.enlarge_window)

        self.progress = QtGui.QProgressBar(self)
        self.progress.setGeometry(0,100,300,10)
        self.progress_btn = QtGui.QPushButton("Increase progressbar", self)
        self.progress_btn.resize(self.progress_btn.sizeHint())
        self.progress_btn.move(0,135)
        self.progress_btn.clicked.connect(self.increase_progress)


        self.show()

    def increase_progress(self):
        self.completed += 1
        self.progress.setValue(self.completed)

    def enlarge_window(self, state):
        if state == QtCore.Qt.Checked:
            self.setGeometry(50,50,600,600)
        else:
            self.setGeometry(50,50, 500, 300)

    def close_application(self):
        choise = QtGui.QMessageBox.question(self,"Extract!", "Are you sure?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

        if choise == QtGui.QMessageBox.Yes:
            print("Closing application...")
            sys.exit()
        else:
            pass

def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())

run()