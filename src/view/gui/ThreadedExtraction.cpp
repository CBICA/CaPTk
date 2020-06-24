#include "ThreadedExtraction.h"

void ThreadedExtraction::run() {
   if (QFile::exists(this->fullPath))
   {
      ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("true").toString());
      ApplicationPreferences::GetInstance()->SerializePreferences();
      ApplicationPreferences::GetInstance()->DisplayPreferences();
      
      QZipReader zr(this->fullPath);
      bool ret = zr.extractAll(this->extractPath);
      
      connect(&zr,SIGNAL(progress(int)),this,SLOT(updateProgressSlot(int)));

      if(ret)
      {
         //after extraction remove the zip
         bool successfullyremoved = QFile::remove(this->fullPath.toStdString().c_str());
      }
   }

   // qDebug() << "Extraction done in background" << this->fullPath << endl;

   emit resultReady(this->appName);
}	

void ThreadedExtraction::updateProgressSlot(int progress) {
   QString msg = "Extracting " + this->appName;
   emit updateProgressSignal(progress, msg.toStdString(), 100;
}