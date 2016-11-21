from __future__ import print_function, division, absolute_import
from stitch_buddy import *
import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *


class EmittingStream(QObject):

	textWritten = pyqtSignal(str)

	def write(self, text):
		self.textWritten.emit(str(text))

class ErrorEmittingStream(QObject):

	textWritten = pyqtSignal(str)

	def write(self, errortext):
		self.textWritten.emit(str(errortext))


class StitchThread(QThread):
	#Runs the stitching in a separate thread in order not to freeze the GUI
	def __init__(self, wellDict, indir, outdir, rescale, RAMstitchFlag):
		QThread.__init__(self)
		self.wellDict = wellDict
		self.indir = indir
		self.outdir = outdir
		self.rescale = rescale
		self.RAMstitchFlag = RAMstitchFlag

	def __del__(self):
		self.quit()

	def run(self):

		if self.RAMstitchFlag:

			stitchWellsInRAM(self.wellDict, self.indir, self.outdir,
							 self.rescale)

		else:
			stitchWellsOnDisk(self.wellDict, self.indir, self.outdir, self.rescale)

class WellRenamer(QWidget):
	"""
	Popup menu for renaming the wells in a wellDict.
	"""

	def __init__(self, wellDict):
		super(WellRenamer, self).__init__()
		self.wellDict = wellDict
		self.inputs = []
		self.initUI()

	def initUI(self):

		fbox = QFormLayout()
		col1_lbl = QLabel("old name:", self)
		col2_lbl = QLabel("new filename:", self)
		fbox.addRow(col1_lbl, col2_lbl)



		for k in self.wellDict.keys():


			label = QLabel(str(k)+":", self)
			txtinput = QLineEdit(str(k))
			fbox.addRow(label, txtinput)
			self.inputs.append((k, label, txtinput))


		qbtn = QPushButton("Cancel", self)
		qbtn.clicked.connect(self.close)
		qbtn.resize(qbtn.sizeHint())

		ok_btn = QPushButton("Ok", self)
		ok_btn.clicked.connect(self.update_welldict)
		fbox.addRow(ok_btn, qbtn)
		self.setLayout(fbox)
		self.setGeometry(500, 300, 200, 300)
		self.setWindowTitle("Rename wells")
		self.setWindowIcon(QIcon("logo.png"))

		self.show()

	def update_welldict(self):
		out = {}
		for entry in self.inputs:
			new_well_name = entry[2].text()
			out[new_well_name] = self.wellDict[entry[0]]
		self.wellDict = out
		self.close()

	def getWelldict(self):
		return self.wellDict




	def __del__(self):
		print("__del__ run")

class AppWindow(QWidget):

	def __init__(self):
		super(AppWindow, self).__init__()

		self.indir = None
		self.outdir = None
		self.wellDict = None
		self.wellDict = filenamesToDict("/Volumes/HDD/Huygens_SYNC/Raw OME files for test/4wells_2x2-mosaik_2-channel_4-frames_test_1/")
		self.rescale = 0.5
		self.RAMstitchFlag = True

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
		lbl1 = QLabel("Load files from this directory:", self)
		btn1_label = "Select input directory"
		self.indir_lbl = QLabel(str(self.indir), self)

		lbl2 = QLabel("Save stitched files to this directory:", self)
		btn2_label = "Select output directory"
		self.outdir_lbl = QLabel(str(self.outdir), self)

		btn3_label = "Rename wells (optional)"

		#Tooltips for the buttons
		btn1_tooltip = "Select the directory where yor OME-TIFF files are " \
					   "stored"
		btn2_tooltip = "Select the directory where your stitched files will " \
					   "be stored"
		btn3_tooltip = "Opens an editor where you can write the names of the" \
					   " stitched wells"

		ramStitchCeckbox_tooltip = "Uncheck ONLY when the stitch operation won't fit in RAM, "\
				"i.e. the size of the files constituting one well plus "\
				"the size of the output is greater than the available RAM"

		#Make buttons work
		btn1 = QPushButton(btn1_label, self)
		btn1.clicked.connect(self.select_indir)
		btn1.setToolTip(btn1_tooltip)
		btn1.resize(btn1.sizeHint())

		btn2 = QPushButton(btn2_label, self)
		btn2.clicked.connect(self.select_outdir)
		btn2.setToolTip(btn1_tooltip)
		btn2.resize(btn2.sizeHint())

		btn3 = QPushButton(btn3_label, self)
		btn3.clicked.connect(self.rename_wells)
		btn3.setToolTip(btn3_tooltip)
		btn3.resize(btn3.sizeHint())

		ramStitchCeckBox = QCheckBox("Stitch in RAM (runs VERY slow if unchecked!)", self)
		ramStitchCeckBox.setToolTip(ramStitchCeckbox_tooltip)
		ramStitchCeckBox.toggle()
		ramStitchCeckBox.stateChanged.connect(self.checkboxState)

		self.rescale_lbl = QLabel("Output will be rescaled with a factor of: "+str(self.rescale))

		rescale_btn = QComboBox(self)
		rescale_btn.addItem("Rescale by 1/2 (default)")
		rescale_btn.addItem("Rescale by 1 (no rescaling)")
		rescale_btn.addItem("Rescale by 1/4")
		rescale_btn.activated[str].connect(self.set_rescale)

		#textbox to log output and errors from component functions
		self.logOutput = QTextEdit()
		self.logOutput.setReadOnly(True)


		qbtn = QPushButton("Quit", self)
		qbtn.clicked.connect(QCoreApplication.instance().quit)
		qbtn.resize(qbtn.sizeHint())

		runButton = QPushButton("Run...", self)
		runButton.clicked.connect(self.run_stitching)
		qbtn.resize(qbtn.sizeHint())

		#Group run and quit buttons into one row
		runQuitBox = QHBoxLayout()
		runQuitBox.addWidget(runButton)
		runQuitBox.addWidget(qbtn)

		#Separator line
		line1 = QFrame()
		line1.setFrameShape(QFrame.HLine)
		line1.setFrameShadow(QFrame.Sunken)


		verticalWidgets = [lbl1, self.indir_lbl, btn1, lbl2, self.outdir_lbl,
						   btn2, line1, btn3, self.rescale_lbl, rescale_btn, ramStitchCeckBox]

		vbox = QVBoxLayout()
		for widget in verticalWidgets:
			vbox.addWidget(widget)

		vbox.addWidget(self.logOutput)
		vbox.addStretch(1)
		vbox.addLayout(runQuitBox)


		self.setLayout(vbox)
		self.setGeometry(300, 300, 900, 500)
		self.setWindowTitle("Stitch buddy")
		self.setWindowIcon(QIcon("logo.png"))
		#self.center()
		self.show()

	def set_rescale(self, text):
		valDict = {"Rescale by 1/2 (default)":0.5, "Rescale by 1 (no rescaling)":None, "Rescale by 1/4":0.25}
		self.rescale = valDict[str(text)]
		self.rescale_lbl.setText("Output will be rescaled with a factor of: "+str(self.rescale))

	def run_stitching(self):
		print("starting")
		try:
			self.wellDict = self.well_rename_popup.getWelldict()
		except:
			pass

		separateThread = StitchThread(self.wellDict, str(self.indir),
									  str(self.outdir), self.rescale,
									  self.RAMstitchFlag)
		separateThread.start()
		print("done!")


	def checkboxState(self, state):

		if state == Qt.Checked:
			self.RAMstitchFlag = True
		else:
			self.RAMstitchFlag = False

		print("RAMstitchFlag set to {}".format(self.RAMstitchFlag))


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
		color = QColor(255, 0, 0)
		self.logOutput.setTextColor(color)

		#Write ouput to log
		self.logOutput.insertPlainText(errortext)

		#Set fontcolor back to black
		color = QColor(0, 0, 0)
		self.logOutput.setTextColor(color)

		#Autoscroll the text
		sb= self.logOutput.verticalScrollBar()
		sb.setValue(sb.maximum())


	def showDialog(self):

		col = QColorDialog.getColor()

		text, ok = QInputDialog.getText(self, "input dialog",
											  "Enter name:")
		if ok:
			self.le.setText(str(text))



	def select_indir(self):
		self.indir = QFileDialog.getExistingDirectory(None,
															'Select input folder...',
															'/Volumes/HDD/Huygens_SYNC/Raw OME files for test/4wells_2x2-mosaik_2-channel_4-frames_test_1/',
															QFileDialog.ShowDirsOnly)

		self.wellDict = filenamesToDict(str(self.indir))
		self.indir_lbl.setText(str(self.indir))
		boldFont = QFont()
		boldFont.setBold(True)
		self.indir_lbl.setFont(boldFont)

	def select_outdir(self):
		self.outdir = QFileDialog.getExistingDirectory(None,
															'Select output folder...',
															'/Volumes/HDD/Huygens_SYNC/Raw OME files for test/testout/',
															QFileDialog.ShowDirsOnly)

		self.outdir_lbl.setText(str(self.outdir))
		boldFont = QFont()
		boldFont.setBold(True)
		self.outdir_lbl.setFont(boldFont)


	def rename_wells(self):

		if self.wellDict == None:
			print("No indir selected!")

		self.well_rename_popup = WellRenamer(self.wellDict)
		self.well_rename_popup.setGeometry(500, 300, 200, 300)
		self.well_rename_popup.setWindowTitle("Rename wells")
		self.well_rename_popup.setWindowIcon(QIcon("logo.png"))

		self.well_rename_popup.show()
		#self.welldict is updated before the stitchingThread is run


	def closeEvent(self, event):

		reply = QMessageBox.question(self, "Message",
										   "Are you sure you want to quit?",
										   QMessageBox.Yes | QMessageBox.No,
										   QMessageBox.Yes)
		if reply == QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()


def main():

	app = QApplication(sys.argv)
	showMe = AppWindow()
	sys.exit(app.exec_())

if __name__ == "__main__":
	main()
