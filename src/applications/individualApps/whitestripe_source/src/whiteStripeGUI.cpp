#include <QApplication>
#include "whiteStripeGUI.h"

//------------------Utility functions----------------------------------
inline bool ifFileExists(QString path) {
  QFileInfo check_file(path);
  // check if file exists and if yes: Is it really a file and no directory?
  if (check_file.exists() && check_file.isFile()) {
    return true;
  }
  else {
    return false;
  }
}
inline bool icasecmp(const string& l, const string& r)
{
  return l.size() == r.size()
    && equal(l.cbegin(), l.cend(), r.cbegin(),
    [](string::value_type l1, string::value_type r1)
  { return toupper(l1) == toupper(r1); });
}
inline bool cmdOptionExists(const vector<string>& allArgs, const string& option)
{
  for (auto &arg : allArgs)
  {
    if (icasecmp(arg, option))
    {
      return true;
    }
  }
  return false;
}
void setStyleSheet()
{
  string styleFileFullPath;

  if (QFile(QApplication::applicationDirPath() + "/captk.qss").exists())
  {
    styleFileFullPath = QApplication::applicationDirPath().toStdString() + "/captk.qss";
  }
  else if (QFile(QApplication::applicationDirPath() + "/../etc/captk.qss").exists())
  {
    styleFileFullPath = QApplication::applicationDirPath().toStdString() + "/../etc/captk.qss";
  }
  else if (QFile(QApplication::applicationDirPath() + "/../Resources/etc/captk.qss").exists())
  {
    styleFileFullPath = QApplication::applicationDirPath().toStdString() + "/../Resources/etc/captk.qss";
  }

  if (!ifFileExists(QString(styleFileFullPath.c_str())))
  {
    styleFileFullPath = "style.qss";
  }
  QFile f(styleFileFullPath.c_str());
	if (!f.exists())
	{
		printf("Unable to set stylesheet, file not found\n");
	}
	else
	{
		f.open(QFile::ReadOnly | QFile::Text);
		QTextStream ts(&f);
		qApp->setStyleSheet(ts.readAll());
	}
}
//-----------------class functions ------------------------------

void WhiteStripeGUI::batchItemDoubleClicked(QListWidgetItem* item)
{
  m_txt_inFile->setText(item->text());
  runAlg();
  cout << "Slot";
}
void WhiteStripeGUI::updateDisplay()
{
  if (m_vtkImg == NULL) return;
  setWindowLevel();
#if VTK_MAJOR_VERSION >= 6
  vtkImageData* display;
#else
  vtkAlgorithmOutput *display;
#endif
  if (m_toggleMask->isChecked() && m_vtkImgMask != NULL)
  {
    displayMask(m_vtkImgMask);
    return;
    display = getDisplayableVtkImg(m_vtkImgMask, true);
  }
  else
  {
    display = getDisplayableVtkImg(m_vtkImg, false);
    m_toggleMask->setChecked(false);
  }
  for (size_t i = 0; i < 3; i++)
  {
    m_views[i]->setEnabled(true);
#if VTK_MAJOR_VERSION >= 6
    m_sliceViewers[i]->SetInputData(display);
#else
    m_sliceViewers[i]->SetInputConnection(display);
#endif
    if (i == 0)
    {
      m_sliceViewers[i]->SetSliceOrientationToXY();
    }
    else if (i == 1)
    {
      m_sliceViewers[i]->SetSliceOrientationToYZ();
    }
    else
    {
      m_sliceViewers[i]->SetSliceOrientationToXZ();
    }
    m_sliceViewers[i]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkInteractorStyleRubberBandZoom::New());
    m_sliceViewers[i]->Render();
    int max = m_sliceViewers[i]->GetSliceMax();
    if (max != m_sliders[i]->maximum())
    {
      m_sliders[i]->setMaximum(max);
      m_sliders[i]->setValue(max / 2);
      m_sliceViewers[i]->SetSlice(max / 2);
    }
    else
    {
      m_sliders[i]->setValue(m_sliders[i]->value());
    }
    m_renderers[i]->ResetCamera();
  }

}
void WhiteStripeGUI::open(QString fileName )
{
  if (fileName.isEmpty())
  {
    QString filter = "MRI Images (*.nii.gz *.nii)";
    fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::currentPath(), filter);
  }
  m_txt_inFile->setText(fileName);
  openMri();
}
void WhiteStripeGUI::addBatchFiles()
{

  QStringList fileNames = QFileDialog::getOpenFileNames(this, tr("Open File"), QDir::currentPath(), "MRI Images (*.nii.gz)");
  m_batchFileList->addItems(fileNames);
  return;
}
void WhiteStripeGUI::save()
{
  if (m_itkImg.IsNull())
  {
    log("Error! No result to save . Please run the algorithm on valid input");
    return;
  }
  QString filters("MRI Images (*.nii.gz");
  QString outFileName = QFileDialog::getSaveFileName(0, "Save file", QDir::currentPath(), filters);

  Utils::writeNifti(m_txt_inFile->text().toStdString(), outFileName.toStdString(), m_itkImg);
  log("Wrote:-" + outFileName.toStdString());
}
void WhiteStripeGUI::imgSliderValuechanged(int val)
{
  QSlider* obj = qobject_cast<QSlider *>(sender());
  size_t index = std::find(m_sliders.begin(), m_sliders.end(), obj) - m_sliders.begin();
  if (index >= 0 && m_sliders.size() > index && m_sliceViewers.size() > index)
  {
    m_sliceViewers[index]->SetSlice(val);
  }
  return;
}
void WhiteStripeGUI::runAlg()
{
  string inFile = m_txt_inFile->text().toStdString();
  string outFile = m_txt_outFile->text().toStdString();

  if (!Utils::fileExist(inFile))
  {
    log("Error. cannot find:" + inFile);
    return;
  }
  else
  {
    log("loading:" + inFile);
    m_itkImg = Utils::readNifti(inFile);
  }

  log("Processing:" + inFile);
  runCore();
  log("Displaying results:" + inFile);
  updateWindowLevel();
  updateDisplay3D();
  log("All done!");


}
void WhiteStripeGUI::runCore()
{
  float twsWidth = m_txt_twsWidth->text().toFloat();
  int sliceStartZ = -1;
  int sliceStopZ = -1;

  if (m_txt_sliceStartZ->isEnabled())
  {
    sliceStartZ = m_txt_sliceStartZ->text().toInt();
    sliceStopZ = m_txt_sliceStopZ->text().toInt();
  }

  int tissuesMax = m_txt_tissuesMax->text().toInt();
  float smoothMax = m_txt_smoothMax->text().toFloat();
  float smoothDelta = m_txt_smoothDelta->text().toFloat();
  int histSize = m_txt_histSize->text().toInt();
  int bT1Int = m_cmb_bT1->currentIndex();
  bool bT1 = (bT1Int == 0);
  m_alg.setParams(twsWidth, sliceStartZ, sliceStopZ, tissuesMax, smoothMax, smoothDelta, histSize, bT1);

  m_itkImg = m_alg.process(m_itkImg, m_itkImg_mask);
  if (m_itkImg.IsNull())
  {
    log("Error processing. Check image validity:");
    m_vtkImg = NULL;
    m_vtkImgMask = NULL;
  }
  else
  {
    m_vtkImg = toVtkImg(m_itkImg);
    m_vtkImgMask = toVtkImg(m_itkImg_mask);
  }

}
void WhiteStripeGUI::runBatch()
{
  int count = m_batchFileList->count();
  if (count == 0)
  {
    log("Error!. no files to process. Please add few!");
    return;
  }
  m_progressBar->setValue(0);
  for (size_t index = 0; index < (uint)count; index++)
  {
    QListWidgetItem * item = m_batchFileList->item(index);
    string inFile = item->text().toStdString();
    string outFile = inFile;
    Utils::replaceSubStr(outFile, ".nii.gz", "");
    outFile = outFile + m_BatchAppendName->text().toStdString() + ".nii.gz";
    m_itkImg = Utils::readNifti(inFile);
    runCore();
    if (m_itkImg.IsNotNull())
    {
      Utils::writeNifti(inFile, outFile, m_itkImg);
      log("Wrote:-" + outFile);
    }
    m_progressBar->setValue((index + 1) * 100 / count);
  }
  m_progressBar->setValue(0);
  log("Done processing batch!");
  return;
}
void WhiteStripeGUI::displayHistChart()
{
  vector<float> mids, origHist, smoothHist;
  vector<int> peakIds;
  int modeId;
  m_alg.getHisInfo(mids, origHist, smoothHist, peakIds, modeId);
  if (mids.empty())
  {
    log("Error no histogram present: First run algorithm to compute histogram!");
    return;
  }

  if (m_hWdg != NULL) delete m_hWdg;
  m_hWdg = new  HistWidget(this);
  const int colCount = 2;
  m_hWdg->setAxis(mids, colCount);
  m_hWdg->addColumn(origHist, "Hist", 1, Scalar(0, 255, 255, 255));
  m_hWdg->addColumn(smoothHist, "Smooth", 2);
  float height = *max_element(smoothHist.begin(), smoothHist.end());
  m_hWdg->plotVirticalLine(mids[modeId], height, "Mode");
  m_hWdg->show();
}
void WhiteStripeGUI::displayAboutAct()
{
  string version = string("WhiteStripe version: ") + WHITESTRIPE_VERSION;
  QMessageBox::information(this, tr(version.c_str()),
    tr("Contact: software@cbica.upenn.edu\n Copyright (c) 2017 University of Pennsylvania. All rights reserved.\n See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html")
    );
}
void WhiteStripeGUI::unloadImage()
{
  m_vtkImg = vtkImageData::New();
  for (size_t i = 0; i < m_sliceViewers.size(); i++) {
    m_sliceViewers[i]->SetInputData(m_vtkImg);
    m_sliceViewers[i]->Render();
  }
}
void WhiteStripeGUI::updateWindowLevel()
{
  double range[2];
  if (m_vtkImg == NULL) return;
  m_vtkImg->GetScalarRange(range);
  m_window = range[1] - range[0];
  m_level = (range[1] + range[0]) * 0.5;
  updateDisplay();
}
void WhiteStripeGUI::toggleMaskClicked()
{
  if (m_toggleMask->isChecked() && m_vtkImgMask == NULL)
  {
    log("Error mask not present: First run algorithm to toggle mask and image display!");
    m_toggleMask->setChecked(false);
  }
  updateDisplay();
}
void WhiteStripeGUI::setWindowLevel()
{
  if (m_vtkImg == NULL) return;
  double window = m_window;
  double level = m_level;
  if (m_toggleMask->isChecked() && m_vtkImgMask != NULL)
  {
    window = 255;
    level = 127.5;
  }
  for (size_t i = 0; i < m_sliceViewers.size(); i++)
  {
    m_sliceViewers[i]->SetColorWindow(window);
    m_sliceViewers[i]->SetColorLevel(level);
    m_views[i]->update();
  }
  return;

}


vtkSmartPointer<vtkImageData> WhiteStripeGUI::toVtkImg(ImageType::Pointer image)
{
  if (image.IsNotNull())
  {
    auto connector = itk::ImageToVTKImageFilter<ImageType>::New();
    connector->SetInput(image);
    connector->Update();
    return connector->GetOutput();
  }
  else
  {
    return vtkSmartPointer<vtkImageData>::New();
  }
}
vtkImageData* WhiteStripeGUI::getDisplayableVtkImg(vtkImageData* vtkImg, bool mask)
{
  if (mask)
  {
    vtkSmartPointer<vtkImageMapToColors> colorMap = vtkSmartPointer<vtkImageMapToColors>::New();
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfColors(2);
    lut->Build();
    lut->SetTableRange(0, 1);
    lut->SetTableValue(0, 0, 0, 0);
    lut->SetTableValue(1, 0, 1, 0);
    
    colorMap->SetLookupTable(lut);
    colorMap->SetInputData(vtkImg);
    //colorMap->SetInputConnection(vtkImg->GetProducerPort());
    colorMap->Update();
    return colorMap->GetOutput();
  }
  return m_vtkImg;
}
void WhiteStripeGUI::displayMask(vtkImageData* image)
{

  vtkSmartPointer<vtkImageMapToColors> colorMap = vtkSmartPointer<vtkImageMapToColors>::New();
  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfColors(2);
  lut->Build();
  lut->SetTableRange(0, 1);
  lut->SetTableValue(0, 0, 0, 0);
  lut->SetTableValue(1, 0, 1, 0);


  colorMap->SetLookupTable(lut);
  colorMap->SetInputData(image);

  for (size_t i = 0; i < 3; i++)
  {
    m_views[i]->setEnabled(true);
    m_sliceViewers[i]->SetInputConnection(colorMap->GetOutputPort());
    if (i == 0)
    {
      m_sliceViewers[i]->SetSliceOrientationToXY();
    }
    else if (i == 1)
    {
      m_sliceViewers[i]->SetSliceOrientationToYZ();
    }
    else
    {
      m_sliceViewers[i]->SetSliceOrientationToXZ();
    }
    m_sliceViewers[i]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkInteractorStyleRubberBandZoom::New());
    m_sliceViewers[i]->Render();
    int max = m_sliceViewers[i]->GetSliceMax();
    if (max != m_sliders[i]->maximum())
    {
      m_sliders[i]->setMaximum(max);
      m_sliders[i]->setValue(max / 2);
      m_sliceViewers[i]->SetSlice(max / 2);
    }
    else
    {
      m_sliders[i]->setValue(m_sliders[i]->value());
    }
    m_renderers[i]->ResetCamera();
  }

}
void WhiteStripeGUI::updateDisplay3D()
{
  vtkActorCollection* actorCollection = m_renderers[3]->GetActors();
  if (actorCollection != NULL)
  {
    actorCollection->InitTraversal();
    for (vtkIdType i = 0; i < actorCollection->GetNumberOfItems(); i++)
    {
      vtkActor* nextActor = actorCollection->GetNextActor();
      m_renderers[3]->RemoveActor(nextActor);
    }
  }

  if (m_vtkImgMask == NULL)
  {
    m_views[3]->update();
    return;
  }
  auto imageDataGeometryFilter = vtkImageDataGeometryFilter::New();
  imageDataGeometryFilter->SetInputData(m_vtkImgMask);
  imageDataGeometryFilter->Update();
  auto thresholdPoints = vtkThresholdPoints::New();
  thresholdPoints->SetInputConnection(imageDataGeometryFilter->GetOutputPort());
  thresholdPoints->ThresholdByUpper(0.5);
  thresholdPoints->Update();

  auto lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfColors(2);
  lut->Build();
  lut->SetTableRange(0, 1);
  lut->SetTableValue(0, 0, 0, 0);
  lut->SetTableValue(1, 0, 1, 0);

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  mapper->SetLookupTable(lut);
  mapper->SetInputConnection(thresholdPoints->GetOutputPort());
  vtkActor* actor = vtkActor::New();
  actor->SetMapper(mapper);


  m_renderers[3]->AddActor(actor);
  m_renderers[3]->ResetCamera();
  m_views[3]->update();
}
void WhiteStripeGUI::openMri()
{
  m_itkImg = Utils::readNifti(m_txt_inFile->text().toStdString());
  m_vtkImg = toVtkImg(m_itkImg);
  m_vtkImgMask = NULL;
  m_itkImg = NULL;
  m_itkImg_mask = NULL;
  updateDisplay();
}
void WhiteStripeGUI::createMenus()
{
  QAction* openAct = new QAction(tr("&Open..."), this);
  openAct->setShortcut(tr("Ctrl+O"));
  connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

  QAction* saveAct = new QAction(tr("&Save Result."), this);
  saveAct->setShortcut(tr("Ctrl+S"));
  connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

  QMenu* fileMenu = new QMenu(tr("&File"), this);
  fileMenu->addAction(openAct);
  fileMenu->addAction(saveAct);

  QAction* displayHistChartAct = new QAction(tr("&Display Histogram..."), this);
  connect(displayHistChartAct, SIGNAL(triggered()), this, SLOT(displayHistChart()));
  QMenu* histMenu = new QMenu(tr("&Histogram"), this);
  histMenu->addAction(displayHistChartAct);

  QAction* displayAboutAct = new QAction(tr("&About WhiteStripe..."), this);
  connect(displayAboutAct, SIGNAL(triggered()), this, SLOT(displayAboutAct()));
  QMenu* aboutMenu = new QMenu(tr("&About"), this);
  aboutMenu->addAction(displayAboutAct);

  menuBar()->addMenu(fileMenu);
  menuBar()->addMenu(histMenu);
  menuBar()->addMenu(aboutMenu);
}

void WhiteStripeGUI::CreateDummyImage(vtkSmartPointer<vtkImageData>& image)
{
  image->SetDimensions(1, 1, 1);
  image->AllocateScalars(VTK_FLOAT, 3);

  int* dims = image->GetDimensions();
  for (int y = 0; y < dims[1]; y++)
  {
    for (int x = 0; x < dims[0]; x++)
    {
      float* pixel = static_cast<float*>(image->GetScalarPointer(x, y, 0));
      pixel[0] = 0;
      pixel[1] = 0;
      pixel[2] = 0;
    }
  }
}
QGridLayout* WhiteStripeGUI::createViews()
{

  m_gridLayout = new QGridLayout();
  for (int i = 0; i < 4; i++)
  {
    QVTKWidget* view = new QVTKWidget();;
    view->setMinimumSize(256, 256);
    view->setEnabled(true);
    if (i < 3)
    {
      vtkImageViewer2* sliceViewer = vtkImageViewer2::New();
      sliceViewer->SetupInteractor(view->GetInteractor());
      sliceViewer->SetRenderWindow(view->GetRenderWindow());
      m_vtkImg = vtkImageData::New();
      CreateDummyImage(m_vtkImg);
      sliceViewer->SetInputData(m_vtkImg);
      m_renderers.push_back(sliceViewer->GetRenderer());
      m_sliceViewers.push_back(sliceViewer);
    }
    else
    {
      vtkRenderer* ren = vtkRenderer::New();
      vtkRenderWindow *renWin = view->GetRenderWindow();
      renWin->AddRenderer(ren);
      m_renderers.push_back(ren);
      view->GetRenderWindow()->GetInteractor()->Initialize();
    }
    m_views.push_back(view);
    QSlider* slider = new QSlider(Qt::Horizontal);
    QObject::connect(slider, SIGNAL(valueChanged(int)), this, SLOT(imgSliderValuechanged(int)));
    m_sliders.push_back(slider);
    view->installEventFilter(this);//Handle double click
  }

  m_gridLayout->addWidget(m_views[0], 0, 0);
  m_gridLayout->addWidget(m_views[1], 0, 1);
  m_gridLayout->addWidget(m_sliders[0], 1, 0);
  m_gridLayout->addWidget(m_sliders[1], 1, 1);
  m_gridLayout->addWidget(m_views[2], 2, 0);
  m_gridLayout->addWidget(m_views[3], 2, 1);
  m_gridLayout->addWidget(m_sliders[2], 3, 0);
  m_gridLayout->addWidget(m_sliders[3], 3, 1);
  return m_gridLayout;
}
QVBoxLayout* WhiteStripeGUI::creatSideBar()
{
  int lblWidth = QLabel().fontMetrics().width("Default Long TextCap");
  QVBoxLayout* sideBarLayout = new QVBoxLayout();
  QScrollArea* algTab = createAlgTab(lblWidth);
  QScrollArea* batchTab = createBatchTab(lblWidth);
  QTabWidget* tabWidget = new QTabWidget();
  tabWidget->setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding));
  tabWidget->addTab(algTab, "SETTINGS");
  //tabWidget->addTab(preProcessTab, "PRE PROCESS");
  tabWidget->addTab(batchTab, "BATCH PROCESSING");
  sideBarLayout->addWidget(tabWidget);


  m_toggleMask = new QCheckBox("Toggle Mask/Image");
  connect(m_toggleMask, SIGNAL(clicked()), this, SLOT(toggleMaskClicked()));
  QPushButton* btnAutoLevel = new QPushButton("Level Display");
  connect(btnAutoLevel, SIGNAL(clicked()), this, SLOT(updateWindowLevel()));

  QPushButton* btnRunWhiteStripe = new QPushButton("Run WhiteStripe");
  connect(btnRunWhiteStripe, SIGNAL(clicked()), this, SLOT(runAlg()));

  QPushButton* btnRunBatch = new QPushButton("Run Batch");
  connect(btnRunBatch, SIGNAL(clicked()), this, SLOT(runBatch()));
  QHBoxLayout* dispLayout = new QHBoxLayout();
  btnAutoLevel->setFixedWidth(lblWidth);
  m_toggleMask->setFixedWidth(lblWidth);
  dispLayout->addWidget(btnAutoLevel);
  dispLayout->addWidget(m_toggleMask);
  sideBarLayout->addLayout(dispLayout);

  QHBoxLayout* controlLayout = new QHBoxLayout();
  btnRunWhiteStripe->setFixedWidth(lblWidth);
  btnRunBatch->setFixedWidth(lblWidth);
  controlLayout->addWidget(btnRunWhiteStripe);
  controlLayout->addWidget(btnRunBatch);
  sideBarLayout->addLayout(controlLayout);

  m_progressLable = new QLabel();
  m_logOutput = new QTextEdit();
  m_progressBar = new QProgressBar();
  QPushButton* logClrButton = new QPushButton("Clear");
  logClrButton->setFixedWidth(QLabel().fontMetrics().width("Clear") + 10);
  connect(logClrButton, SIGNAL(clicked()), m_progressLable, SLOT(clear()));
  connect(logClrButton, SIGNAL(clicked()), m_logOutput, SLOT(clear()));
  connect(logClrButton, SIGNAL(clicked()), m_progressBar, SLOT(reset()));
  m_progressBar->setFixedWidth(2 * lblWidth - logClrButton->width());
  QHBoxLayout* hlYoutLogCap = new QHBoxLayout();
  hlYoutLogCap->addWidget(m_progressBar);
  hlYoutLogCap->addWidget(logClrButton);
  sideBarLayout->addWidget(m_progressLable);
  sideBarLayout->addLayout(hlYoutLogCap);

  m_logOutput->setReadOnly(true);
  m_logOutput->setMinimumWidth(2 * lblWidth); //Horizondal Scrolling is not working here TBD 
  m_logOutput->setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Expanding));
  //cout << "Width= " << m_logOutput->width() << " " << 2 * lblWidth;
  m_logOutput->setLineWrapMode(QTextEdit::NoWrap);
  m_logOutput->setHorizontalScrollBarPolicy(Qt::ScrollBarPolicy::ScrollBarAlwaysOn);
  sideBarLayout->addWidget(m_logOutput);
  return sideBarLayout;
}
QScrollArea* WhiteStripeGUI::createBatchTab(const int lblWidth)
{
  QScrollArea* scrollArea = new QScrollArea();
  QWidget* scrollWdg = new QWidget();
  QVBoxLayout*  barLayout = new QVBoxLayout(scrollWdg);

  m_batchFileList = new QListWidget();
  m_batchFileList->setMaximumWidth(2 * lblWidth);
  connect(m_batchFileList, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(batchItemDoubleClicked(QListWidgetItem*)));

  QHBoxLayout* controlsLayout = new QHBoxLayout();
  QHBoxLayout* appendLayout = new QHBoxLayout();
  QLabel* appendLabel = new QLabel("Append result:");
  m_BatchAppendName = new QLineEdit("_norm");
  QPushButton* btnAdd = new QPushButton("Add");
  connect(btnAdd, SIGNAL(clicked()), this, SLOT(addBatchFiles()));
  QPushButton* btnClear = new QPushButton("Clear");
  connect(btnClear, SIGNAL(clicked()), m_batchFileList, SLOT(clear()));

  appendLabel->setFixedWidth(lblWidth);
  m_BatchAppendName->setFixedWidth(lblWidth);
  btnAdd->setFixedWidth(lblWidth);
  btnClear->setFixedWidth(lblWidth);

  controlsLayout->addWidget(btnAdd);
  controlsLayout->addWidget(btnClear);
  appendLayout->addWidget(appendLabel);
  appendLayout->addWidget(m_BatchAppendName);
  barLayout->addLayout(appendLayout);
  barLayout->addLayout(controlsLayout);
  barLayout->addWidget(m_batchFileList);
  scrollArea->setWidget(scrollWdg);
  return scrollArea;
}
QScrollArea* WhiteStripeGUI::createAlgTab(const int lblWidth)
{

  QScrollArea* scrollArea = new QScrollArea();
  QWidget* scrollWdg = new QWidget();
  QVBoxLayout*  barLayout = new QVBoxLayout(scrollWdg);
  QHBoxLayout* settingsLayout0_1 = new QHBoxLayout();
  m_txt_inFile = new QLineEdit("in.nii.gz");
  QPushButton* btnBrowse = new QPushButton("Input File:");
  connect(btnBrowse, SIGNAL(clicked()), this, SLOT(open()));

  btnBrowse->setFixedWidth(lblWidth);
  m_txt_inFile->setFixedWidth(lblWidth);
  settingsLayout0_1->addWidget(btnBrowse);
  settingsLayout0_1->addWidget(m_txt_inFile);


  QHBoxLayout* settingsLayout0_2 = new QHBoxLayout();
  m_txt_outFile = new QLineEdit("out.nii.gz");
  QLabel* label = new QLabel("Output File:");
  label->setFixedWidth(lblWidth);
  m_txt_outFile->setFixedWidth(lblWidth);
  settingsLayout0_2->addWidget(label);
  settingsLayout0_2->addWidget(m_txt_outFile);

  QHBoxLayout* settingsLayout1 = new QHBoxLayout();
  m_cmb_bT1 = new QComboBox();
  m_cmb_bT1->addItem("T1");
  m_cmb_bT1->addItem("T2");
  label = new QLabel("Mri Type:");
  label->setFixedWidth(lblWidth);
  m_cmb_bT1->setFixedWidth(lblWidth);
  settingsLayout1->addWidget(label);
  settingsLayout1->addWidget(m_cmb_bT1);

  QHBoxLayout* settingsLayout2 = new QHBoxLayout();
  m_txt_histSize = new QLineEdit("2000");
  m_txt_histSize->setValidator(new QIntValidator(100, 5000, this));
  label = new QLabel("HistBin Size:");
  label->setFixedWidth(lblWidth);
  m_txt_histSize->setFixedWidth(lblWidth);
  settingsLayout2->addWidget(label);
  settingsLayout2->addWidget(m_txt_histSize);


  QHBoxLayout* settingsLayout3 = new QHBoxLayout();
  m_txt_twsWidth = new QLineEdit("0.05");
  m_txt_twsWidth->setValidator(new QDoubleValidator(0.0, 1.0, 2, this));
  label = new QLabel("Wstripe Radius:");
  label->setFixedWidth(lblWidth);
  m_txt_twsWidth->setFixedWidth(lblWidth);
  settingsLayout3->addWidget(label);
  settingsLayout3->addWidget(m_txt_twsWidth);



  QHBoxLayout* settingsLayout4 = new QHBoxLayout();
  m_txt_sliceStartZ = new QLineEdit("80");
  m_txt_sliceStartZ->setValidator(new QIntValidator(0, 500, this));
  label = new QLabel("Slice Start(Z)");
  label->setFixedWidth(lblWidth);
  m_txt_sliceStartZ->setFixedWidth(lblWidth);
  settingsLayout4->addWidget(label);
  settingsLayout4->addWidget(m_txt_sliceStartZ);

  QHBoxLayout* settingsLayout4_2 = new QHBoxLayout();
  m_txt_sliceStopZ = new QLineEdit("120");
  m_txt_sliceStopZ->setValidator(new QIntValidator(0, 500, this));
  label = new QLabel("Slice Stop(Z)");
  label->setFixedWidth(lblWidth);
  m_txt_sliceStopZ->setFixedWidth(lblWidth);
  settingsLayout4_2->addWidget(label);
  settingsLayout4_2->addWidget(m_txt_sliceStopZ);
  QVBoxLayout *vbox = new QVBoxLayout;
  vbox->addLayout(settingsLayout4);
  vbox->addLayout(settingsLayout4_2);
  QGroupBox *axialSlicing = new QGroupBox("Enable Axial Slicing?");
  axialSlicing->setCheckable(true);
  axialSlicing->setChecked(true);
  axialSlicing->setLayout(vbox);

  QHBoxLayout* settingsLayout5 = new QHBoxLayout();
  m_txt_tissuesMax = new QLineEdit("5");
  m_txt_tissuesMax->setValidator(new QIntValidator(2, 20, this));
  label = new QLabel("Max Tissues:");
  label->setFixedWidth(lblWidth);
  m_txt_tissuesMax->setFixedWidth(lblWidth);
  settingsLayout5->addWidget(label);
  settingsLayout5->addWidget(m_txt_tissuesMax);

  QHBoxLayout* settingsLayout6 = new QHBoxLayout();
  m_txt_smoothMax = new QLineEdit("10.0");
  m_txt_smoothMax->setValidator(new QDoubleValidator(1.0, 30.0, 2, this));
  label = new QLabel("Max Smoothing:");
  label->setFixedWidth(lblWidth);
  m_txt_smoothMax->setFixedWidth(lblWidth);
  settingsLayout6->addWidget(label);
  settingsLayout6->addWidget(m_txt_smoothMax);

  QHBoxLayout* settingsLayout7 = new QHBoxLayout();
  m_txt_smoothDelta = new QLineEdit("0.50");
  m_txt_smoothDelta->setValidator(new QDoubleValidator(0.1, 2.0, 2, this));
  label = new QLabel("Smoothing Delta:");
  label->setFixedWidth(lblWidth);
  m_txt_smoothDelta->setFixedWidth(lblWidth);
  settingsLayout7->addWidget(label);
  settingsLayout7->addWidget(m_txt_smoothDelta);

  barLayout->addLayout(settingsLayout0_1);
  barLayout->addLayout(settingsLayout0_2);
  barLayout->addLayout(settingsLayout1);
  barLayout->addLayout(settingsLayout2);
  barLayout->addLayout(settingsLayout3);
  barLayout->addLayout(settingsLayout5);
  barLayout->addLayout(settingsLayout6);
  barLayout->addLayout(settingsLayout7);
  barLayout->addWidget(axialSlicing);
  scrollArea->setWidget(scrollWdg);
  return scrollArea;
}

float  WhiteStripeGUI::parseProgress(string original)
{
  float  progress = -1.0;
  size_t s = original.find("Progress=");
  if (s == string::npos) return progress;
  size_t e = original.find("%", s);
  if (e == string::npos) return progress;
  string sub = original.substr(s + 1, e - s - 1);
  progress = atof(sub.c_str());
  return progress;
}
void WhiteStripeGUI::log(string msg, bool bold, string colour)
{
  float progress = parseProgress(msg);
  if (progress > 0.0f)
  {
    m_progressBar->setValue(int(progress));
    QCoreApplication::processEvents();
    return;
  }

  if (!colour.empty())
  {
    msg = "<font color=\"" + colour + "\"> " + msg + " </font>";
  }
  else
  {
    Utils::replaceSubStr(msg, "Finished", " <font color = \"Green\"> Finished </font>");
    Utils::replaceSubStr(msg, "Success", "<font color=\"Green\"> Success </font>");
    Utils::replaceSubStr(msg, "Error", "<font color=\"Red\"> Error </font>");
    Utils::replaceSubStr(msg, "Cannot", "<font color=\"Yellow\"> Cannot </font>");
    Utils::replaceSubStr(msg, "CORE:", "<font color=\"Maroon\"> CORE: </font>");
  }
  msg + msg + "<br/>";
  msg = QTime::currentTime().toString().toStdString() + " : " + msg;
  if (bold)
  {
    msg = "<b>" + msg + "</b>";
  }
  m_logOutput->textCursor().insertHtml(QString::fromStdString((msg)));
  m_logOutput->textCursor().insertText("\n");
  m_logOutput->ensureCursorVisible();
  QCoreApplication::processEvents();
}
void WhiteStripeGUI::maxDblClick(size_t index)
{
  bool doMin = m_views[0]->isHidden() || m_views[1]->isHidden();//Check atleaset one of them is hidden
  for (size_t i = 0; i < m_views.size(); i++)
  {
    if (doMin)
    {
      m_views[i]->show();
      m_sliders[i]->show();
    }
    else if (i != index)
    {
      m_views[i]->hide();
      m_sliders[i]->hide();
    }
  }
  m_gridLayout->update();
}
bool WhiteStripeGUI::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::MouseButtonDblClick)
  {
    size_t index = std::find(m_views.begin(), m_views.end(), obj) - m_views.begin();
    if (index >= 0 && index< m_views.size())
    {
      maxDblClick(index);
    }
    event->ignore();
    //cout << "\n" << "Mouse Button Double Clicked";
    return true;
  }
  return false;
}
//--------------------------------------------

int main(int argc, char *argv[]) {
  
  std::vector<string> allArgs(argv, argv + argc);
  if (cmdOptionExists(allArgs, "-v") || cmdOptionExists(allArgs, "-version"))
  {
    cout << "WhiteStripe version: " << WHITESTRIPE_VERSION << endl;
    return EXIT_SUCCESS;
  }
  if (cmdOptionExists(allArgs, "-h") || cmdOptionExists(allArgs, "-help"))
  {
    cout << "This is a GUI application. Command Line arguments not supported!" << WHITESTRIPE_VERSION << endl;
    return EXIT_SUCCESS;
  }
  if (cmdOptionExists(allArgs, "-test"))
  {
    cout << "Dummy test for checking build sanity"<< endl;
    return EXIT_SUCCESS;
  }

  QApplication app(argc, argv);
	setStyleSheet();
	WhiteStripeGUI *window = new WhiteStripeGUI(NULL);
	window->show();
  return app.exec();
}
