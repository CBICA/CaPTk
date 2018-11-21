#! /usr/bin/env python
import os, sys, time
try:
  from PyQt4 import QtCore, QtGui
  import numpy as np
  from datetime import datetime
  import vtk
  from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
  sys.path.append(os.path.abspath(os.path.join('..', 'libexec'))) #import ConfettiCore in linux
  import ConfettiCore
  import vtk.util.numpy_support as vnps
  import struct
except:
  print "Failed to import required modules"
  sys.exit(1)


class DtiUtils:
    def numpy_to_vtk_points(self, points):
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(vtk.util.numpy_support.numpy_to_vtk(np.asarray(points), deep=True))
        return vtk_points
    def numpy_to_vtk_colors(self, colors):
        vtk_colors = vtk.util.numpy_support.numpy_to_vtk(np.asarray(colors), deep=True,
                                     array_type=vtk.VTK_UNSIGNED_CHAR)
        return vtk_colors
    def lines_to_vtk_polydata(self,lines):
        points_array = np.vstack(lines)
        nb_lines = len(lines)
        nb_points = len(points_array)
        lines_range = range(nb_lines)
        lines_array = []
        points_per_line = np.zeros([nb_lines], np.intp)
        current_position = 0
        for i in lines_range:
            current_len = len(lines[i])
            points_per_line[i] = current_len

            end_position = current_position + current_len
            lines_array += [current_len]
            lines_array += range(current_position, end_position)
            current_position = end_position
        lines_array = np.array(lines_array)
        vtk_points = self.numpy_to_vtk_points(points_array)

        vtk_lines = vtk.vtkCellArray()
        vtk_lines.GetData().DeepCopy(vtk.util.numpy_support.numpy_to_vtk(lines_array))
        vtk_lines.SetNumberOfCells(nb_lines)


        cols_arr = self.line_colors(lines)
        colors_mapper = np.repeat(lines_range, points_per_line, axis=0)
        vtk_colors = self.numpy_to_vtk_colors(255 * cols_arr[colors_mapper])
        vtk_colors.SetName("Colors")
        # Create the poly_data
        poly_data = vtk.vtkPolyData()
        poly_data.SetPoints(vtk_points)
        poly_data.SetLines(vtk_lines)
        poly_data.GetPointData().SetScalars(vtk_colors)
        return poly_data
    def get_streamlines(self,line_polydata):
        lines_vertices = vtk.util.numpy_support.vtk_to_numpy(line_polydata.GetPoints().GetData())
        lines_idx = vtk.util.numpy_support.vtk_to_numpy(line_polydata.GetLines().GetData())
    
        lines = []
        current_idx = 0
        while current_idx < len(lines_idx):
            line_len = lines_idx[current_idx]
            #print line_len
            next_idx = current_idx + line_len + 1 
            line_range = lines_idx[current_idx + 1: next_idx]
            #print line_range
            lines += [lines_vertices[line_range]]
            current_idx = next_idx
        return lines
    def get_points(self,point_polydata):
        vertices = vtk.util.numpy_support.vtk_to_numpy(point_polydata.GetPoints().GetData())
        return vertices
    def set_input(self, vtk_object, inp):
        if isinstance(inp, vtk.vtkPolyData) \
           or isinstance(inp, vtk.vtkImageData):
            if vtk.VTK_MAJOR_VERSION <= 5:
                vtk_object.SetInput(inp)
            else:
                vtk_object.SetInputData(inp)
        elif isinstance(inp, vtk.vtkAlgorithmOutput):
            vtk_object.SetInputConnection(inp)

        vtk_object.Update()
        return vtk_object
    def line_colors(self,streamlines):
        col_list =[]
        for  streamline in streamlines:
            orient= streamline[-1] - streamline[0]            
            orient= np.abs(orient/np.sqrt(np.sum(np.square(orient))))
            col_list.append(orient)
        return np.vstack(col_list)

    def getStreamLines(self,lines, colors=None, opacity=1, linewidth=1):
        poly_data = self.lines_to_vtk_polydata(lines)
        poly_mapper = self.set_input(vtk.vtkPolyDataMapper(), poly_data)
        poly_mapper.ScalarVisibilityOn()
        poly_mapper.SetScalarModeToUsePointFieldData()
        poly_mapper.SelectColorArray("Colors")
        poly_mapper.Update()

        actor = vtk.vtkActor()
        actor.SetMapper(poly_mapper)
        actor.GetProperty().SetLineWidth(linewidth)
        actor.GetProperty().SetOpacity(opacity)
        return actor
    def testDtiUtils(self):
        lines = [np.array([[-1, 0, 0.], [1, 0, 0.]]),
                 np.array([[-1, 1, 0.], [1, 1, 0.]])]
        steamActor = self.getStreamLines(lines)
        return steamActor

    def readFibersFromFile(self,bfloatFile):
        f = open(bfloatFile, "rb")
        lines =[]
        while True:
            byte = f.read(4)
            if len(byte)<4:
                break
            numPoints = int(struct.unpack('>f', byte)[0])
            byte = f.read(4)
            seedPoint = int(struct.unpack('>f', byte)[0])
            coords = np.zeros([numPoints, 3])
            for i in range(0, numPoints):
                byte = f.read(4)
                x = float(struct.unpack('>f', byte)[0])
                byte = f.read(4)
                y = float(struct.unpack('>f', byte)[0])
                byte = f.read(4)
                z = float(struct.unpack('>f', byte)[0])
                coords[i, :] = [x,y,z]
            lines.append(np.array(coords, dtype=np.float32))
        f.close()
        return np.array(lines)

class DTIApp(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(DTIApp, self).__init__(parent)
        self._isValid= False

        self.api = ConfettiCore.ConfettiApi()
        self.api.registerCppCallback()

        #Create Main GUI
        self.frame = QtGui.QFrame()
        self.mainLayout = QtGui.QHBoxLayout()
        self.createViews()
        self.createMenubar()
        self.createSideBar()
        self.mainLayout.addLayout(self.sideBarLayout)
        self.setCentralWidget(self.frame)
        self.frame.setLayout(self.mainLayout)
        
        # Other
        self.setMyStyleSheet()
        self.setWindowTitle("DTI ConfettiViewer")
        self.log("Log active:",True,"Blue")
        self.api.setListner(self.log)
        self.dtiActor =None
        if self.api.CoreModule is None:
            self.log("Failed to load Core Module ")
        self._initialized = False
        self.show()
        
        

    def showEvent(self, evt):
        if not self._initialized:
            self.iren.Initialize()
            self.startTimer(30)
            self._initialized = True

    def setMyStyleSheet(self):
        if getattr(sys, 'frozen', False):
            application_path = os.path.dirname(sys.executable)
        elif __file__:
            application_path = os.path.dirname(__file__)
        baseDir=os.path.join(application_path, '..' )
        styleFileFullPath =""
        for root, dirs, files in os.walk(baseDir):
            for file in files:
                if "captk.qss" in file:
                     styleFileFullPath= os.path.join(root, file)
        if os.path.isfile(styleFileFullPath) :
            with open(styleFileFullPath, 'r') as myfile:
                data=myfile.read()
                QtGui.qApp.setStyleSheet(data)
        else:
          self.setStyleSheet("""
                    *
                    {
                        background-color: rgb(40, 40, 40);
                        alternate-background-color: rgb(25, 25, 25);
                        color: rgb(100 ,100, 100);
                        selection-background-color: rgb(200, 200, 200);
                    }
                    QMenuBar::item
                    {
                        background-color: rgb(60, 60, 60);
                    }
                    """)

    def createViews(self):
        self.view= QVTKRenderWindowInteractor(self.frame)
        self.view.setMinimumSize(512,512)
        self.renderer = vtk.vtkRenderer()
        self.view.GetRenderWindow().AddRenderer(self.renderer)
        self.view.update()
        self.iren = self.view.GetRenderWindow().GetInteractor()
        #self.view.GetRenderWindow().GetInteractor().Initialize()
        self.mainLayout.addWidget(self.view)


    def open(self):
        dialog = QtGui.QFileDialog(self)
        dialog.setWindowTitle('Open DTI Track file')
        dialog.setNameFilter('Images (*.trk *.Bfloat )')
        dialog.setFileMode(QtGui.QFileDialog.ExistingFile)
        if dialog.exec_() == QtGui.QDialog.Accepted:
            filename = dialog.selectedFiles()[0]
            self.openDtiTrack(str(filename));
    def close(self):
        if self.dtiActor :
            self.renderer.RemoveActor(self.dtiActor)
            self.renderer.ResetCamera()
            self.view.update()
    def closeEvent(self,event):
        self.close()
        self.deleteLater()


    def browseTdiFile(self):
        file = QtGui.QFileDialog.getOpenFileName(self, "Open TDI File","","CSV (*.csv)")
        if (file):
          self.txtTDI_File.setText(str(file))

    def browseOutDir(self):
        dir = QtGui.QFileDialog.getExistingDirectory(self, "Select Directory")
        if (dir):
          self.outputDir.setText(str(dir))
    def browseStreamLineFile(self):
        fileName = QtGui.QFileDialog.getOpenFileName(self,
                "Open Streamline", '', "Streamline Files (*.Bfloat)")
        if fileName:
          self.txtStreamLineFile.setText(str(fileName))

    def creatConfettiTab(self, lblWidth):
        scrollArea = QtGui.QScrollArea()
        scrollWdg= QtGui.QWidget()
        barLayout = QtGui.QVBoxLayout(scrollWdg)


        streamLineLayout = QtGui.QHBoxLayout()
        self.txtStreamLineFile = QtGui.QLineEdit("")
        btnStreamlineFile =QtGui.QPushButton("Streamline File")
        btnStreamlineFile.clicked.connect(self.browseStreamLineFile)
        btnStreamlineFile.setFixedWidth(lblWidth)
        self.txtStreamLineFile.setFixedWidth(lblWidth)
        streamLineLayout.addWidget(btnStreamlineFile)
        streamLineLayout.addWidget(self.txtStreamLineFile)

        tdiLayout = QtGui.QHBoxLayout()
        self.txtTDI_File = QtGui.QLineEdit("")
        btnDataFile =QtGui.QPushButton("TDI File")
        btnDataFile.setFixedWidth(lblWidth)
        self.txtTDI_File.setFixedWidth(lblWidth)
        tdiLayout.addWidget(btnDataFile)
        tdiLayout.addWidget(self.txtTDI_File)
        btnDataFile.clicked.connect(self.browseTdiFile)

        outputLayout = QtGui.QHBoxLayout()
        self.outputDir = QtGui.QLineEdit("")
        btnOutDir =QtGui.QPushButton("Output Directory")
        btnOutDir.setFixedWidth(lblWidth)
        self.outputDir.setFixedWidth(lblWidth)
        outputLayout.addWidget(btnOutDir)
        outputLayout.addWidget(self.outputDir)
        btnOutDir.clicked.connect(self.browseOutDir)

       
        barLayout.addLayout(streamLineLayout)
        barLayout.addLayout(tdiLayout)
        barLayout.addLayout(outputLayout)

        outputLayout = self.creatReviewLayout(lblWidth)
        barLayout.addLayout(outputLayout)
        scrollArea.setWidget(scrollWdg)
        return scrollArea

    def listViewClicked(self, item):
        path = str(item.data(QtCore.Qt.UserRole).toPyObject())
        if path.endswith(".Bfloat") or path.endswith(".trk") :
            self.openDtiTrack(str(path))

    def populateListView(self):
        self.listWidget.clear()
        projectDir =str(self.outputDir.text())
        for dirpath, dirnames, filenames in os.walk(projectDir):
            for filename in [f for f in filenames if f.endswith(".trk") or f.endswith(".Bfloat")]:
                if "tract_" in filename:
                    path= os.path.join(dirpath, filename)
                    item =QtGui.QListWidgetItem()
                    item.setData(QtCore.Qt.UserRole, str(path));
                    item.setText(filename);
                    item.setToolTip(path)
                    self.listWidget.addItem(item)

    def creatReviewLayout(self, width):
        barLayout = QtGui.QVBoxLayout()
        self.listWidget = QtGui.QListWidget()
        self.listWidget.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding);
        self.listWidget.setMaximumWidth(2*width)
        self.listWidget.itemDoubleClicked.connect(self.listViewClicked)
        refreshBtn= QtGui.QPushButton("Refresh")
        refreshBtn.clicked.connect(self.populateListView)
        headerLayout = QtGui.QHBoxLayout()
        headerLayout.addWidget(QtGui.QLabel("Output"))
        headerLayout.addStretch()
        headerLayout.addWidget(refreshBtn)

        barLayout.addLayout(headerLayout)
        barLayout.addWidget(self.listWidget)
        return barLayout

    def createSideBar(self):
            lblWidth = QtGui.QLabel().fontMetrics().width("Default Long TextCap.....");
            self.sideBarLayout = QtGui.QVBoxLayout()
            confettiTab = self.creatConfettiTab(lblWidth)
            tabWidget = QtGui.QTabWidget(self)
            tabWidget.setFixedWidth(2*lblWidth+47)
            tabWidget.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding);
            tabWidget.addTab(confettiTab, "CONFETTI SETTINGS")
            self.sideBarLayout.addWidget(tabWidget)
            self.progressLable= QtGui.QLabel(":")
            self.logOutput = QtGui.QTextEdit(self)

            self.progressBar= QtGui.QProgressBar()
            logClrButton= QtGui.QPushButton("Clear Log")
            logClrButton.setFixedWidth(QtGui.QLabel().fontMetrics().width("  Clear Log "));
            logClrButton.clicked.connect(self.progressLable.clear)
            logClrButton.clicked.connect(self.logOutput.clear)
            logClrButton.clicked.connect(self.progressBar.reset)
            self.progressBar.setFixedWidth(2*lblWidth-logClrButton.width())
            
            self.btnRunAll= QtGui.QPushButton("Run Confetti")
            self.btnRunAll.clicked.connect(self.runAll)
            self.btnRunAll.setSizePolicy(QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed)
            btnLayout = QtGui.QHBoxLayout()
            btnLayout.addWidget(self.progressLable)
            btnLayout.addWidget(self.btnRunAll)
            self.sideBarLayout.addLayout(btnLayout)

            wd= tabWidget.width()
            hlYoutLogCap = QtGui.QHBoxLayout()
            hlYoutLogCap.addWidget(self.progressBar)
            hlYoutLogCap.addWidget(logClrButton)
            self.sideBarLayout.addLayout(hlYoutLogCap)
            
            self.logOutput.setReadOnly(True)
            self.logOutput.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
            self.logOutput.setLineWrapMode(QtGui.QTextEdit.NoWrap)
            self.logOutput.setFixedWidth(2*lblWidth+50)

            self.sideBarLayout.addWidget(self.logOutput)

    def createMenubar(self):
        # File Menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.openAct = QtGui.QAction("&Open...", self,shortcut=QtGui.QKeySequence.Open,statusTip="Open an existing file", triggered=self.open)
        self.closeAct = QtGui.QAction("&Close", self,shortcut=QtGui.QKeySequence.Close, statusTip="Close project", triggered=self.close)
        self.menuFile.addAction(self.openAct)
        self.menuFile.addAction(self.closeAct)
        self.menuHelp= self.menuBar().addMenu("&Help")
        self.helpAct = QtGui.QAction("&About", self,shortcut=QtGui.QKeySequence.Open,statusTip="About information", triggered=self.aboutDialog)
        self.menuHelp.addAction(self.helpAct)

    def notImplementedDialog(self):
        QtGui.QMessageBox.information(self, "Not Implemented", "This functionality is not yet implemented!")
    def aboutDialog(self):
        QtGui.QMessageBox.information(self, "About", "ConfettiViewer  Beta 1.2")
    
    def openDtiTrack(self,fileName):
        self.setCursor(QtCore.Qt.WaitCursor)
        actors= self.renderer.GetActors()
        if self.dtiActor :
            self.renderer.RemoveActor(self.dtiActor)

        streamlines =DtiUtils().readFibersFromFile(fileName)

        self.dtiActor = DtiUtils().getStreamLines(streamlines)
        self.renderer.AddActor(self.dtiActor)
        self.renderer.ResetCamera()
        self.view.update()
        self.setCursor(QtCore.Qt.ArrowCursor)
    def showWarning(self,title,message):
        QtGui.QMessageBox.warning(self,title,message,
                    QtGui.QMessageBox.Cancel, QtGui.QMessageBox.NoButton,
                    QtGui.QMessageBox.NoButton)

    def findBundleAtlasDir(self):
        application_path=""
        if getattr(sys, 'frozen', False):
            application_path = os.path.dirname(sys.executable)
        elif __file__:
            application_path = os.path.dirname(__file__)
        return os.path.join(application_path+"/ConfettiBundleAtlas/")

    def runGenConnectivity(self):
        bundleAtlasDir =self.findBundleAtlasDir()
        if not os.path.isdir(bundleAtlasDir):
            self.showWarning("Error", "Cannot find Bundle Atlas Directory at :"+bundleAtlasDir)
            return False
        self.progressLable.setText("Generating connectivity..")
        tdiFile = str(self.txtTDI_File.text())
        if not os.path.isfile(tdiFile):
            self.showWarning("Error", "Cannot find TDI file at :"+tdiFile)
            return False
        outDir =str(self.outputDir.text())
        if not os.path.isdir(outDir):
            self.showWarning("Error", "Cannot find output Directory at :"+outDir)
            return False
        signatureFileName =outDir+ "/signatures.csv"
        fiberFile =str(self.txtStreamLineFile.text())
        self.api.genConnectivitySig(tdiFile, signatureFileName,fiberFile)
        self.progressBar.reset()
        self.progressLable.setText("Finished generating connectivity.")
        return True

    def runGenCluster(self):
        self.progressLable.setText("Generating Clusters..")
        outDir =self.outputDir.text()
        signatureFileName=str(outDir+ "/signatures.csv")
        tempRoot =self.findBundleAtlasDir()
        outFile = str(outDir+ "/clusterIDs.csv")
        self.api.genAdaptiveCluster(signatureFileName, tempRoot, outFile)
        self.progressBar.reset()
        self.progressLable.setText("Finished generating clusters.")
        return True

    def runExtractTracts(self):
        self.progressLable.setText("Extracting Tracts..")
        projectDir =self.outputDir.text()
        outDIR =str(projectDir+ "/tract_")
        tempRoot = self.findBundleAtlasDir()
        inputIDFile = tempRoot + "/tract_ids.csv";
        inputNamesFile = tempRoot + "/tract_names.csv";
        inputLabelFile = str(projectDir+ "/clusterIDs.csv")
        if not os.path.isfile(inputIDFile) or not os.path.isfile(inputNamesFile)   :
            self.log("Error reading:"+inputIDFile)
            self.log("Error reading:"+inputNamesFile)
            self.showWarning("Error", "Contents missing at at Bundle Atlas Directory:"+tempRoot)
            return False
        if  not os.path.isfile(inputLabelFile):
            self.log("Error reading:"+inputLabelFile)
            self.showWarning("Error", "Cannoot read :"+inputLabelFile)
            return False
            
        inputBfloatFile = str(self.txtStreamLineFile.text())
        if not os.path.isfile(inputBfloatFile):
            self.showWarning("Error", "Cannot find Streamline file at:"+inputBfloatFile)
            return False

        self.api.extractTracts(outDIR,tempRoot,inputIDFile,inputNamesFile,inputLabelFile, inputBfloatFile)
        self.progressBar.reset()
        self.progressLable.setText("Finished Extracting Tracts.")
        return True

    def runAll(self):
        res= self.runGenConnectivity()
        if res:
            res= self.runGenCluster()
        if res:
            res= self.runExtractTracts()
        self.populateListView()
        self.progressBar.reset()
        self.progressLable.setText("Finished running full pipeline")

    def parseProgress(self, s):
        try:
            start = s.index( "Progress=" ) + len( "Progress=" )
            end = s.index( "%", start )
            return float(s[start:end])
        except ValueError:
            return None

    def log(self, msg, bold=False, colour =None):
          progress= self.parseProgress(msg)
          if progress:
             self.progressBar.setValue(int(progress))
             QtCore.QCoreApplication.processEvents()
             return 

          if colour is not None:
            msg = "<font color=\"%s\"> "% colour + msg + " </font>"
          msg =msg.replace("Finished", "<font color=\"Green\"> Finished </font>")
          msg =msg.replace("Success", "<font color=\"Green\"> Success </font>")
          msg =msg.replace("Error", "<font color=\"Red\"> Error </font>")
          msg =msg.replace("Cannot", "<font color=\"Yellow\"> Cannot </font>")
          msg =msg.replace("CORE:", "<font color=\"Maroon\"> CORE: </font>")
          timeStr = datetime.now().strftime('%H:%M:%S') + " : " + msg
          if bold:
            timeStr = "<b>" + timeStr+ "<br> </b>"
            self.logOutput.textCursor().insertHtml(timeStr)
          else :
            timeStr = timeStr+ "<br>" 
            self.logOutput.textCursor().insertHtml(timeStr)
          self.logOutput.ensureCursorVisible()
          QtCore.QCoreApplication.processEvents()

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    player = DTIApp()
    player.show()
    sys.exit(app.exec_())
