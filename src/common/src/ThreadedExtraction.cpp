#include "ThreadedExtraction.h"

void ThreadedExtraction::run() {
   if (QFile::exists(this->fullPath))
   {
      ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("true").toString());
      ApplicationPreferences::GetInstance()->SerializePreferences();
      // ApplicationPreferences::GetInstance()->DisplayPreferences();
      
      QZipReader zr(this->fullPath);
      connect(&zr,SIGNAL(progress(int)),this,SLOT(updateProgressSlot(int)));

      bool ret = zr.extractAll(this->extractPath);
      
      if(ret)
      {
         //after extraction remove the zip
         bool successfullyremoved = QFile::remove(this->fullPath.toStdString().c_str());
         emit resultReady(ret);
      }
   }

   // qDebug() << "Extraction done in background" << this->fullPath << endl;


}	

void ThreadedExtraction::updateProgressSlot(int progress) {
   
   QString msg = "Extracting " + this->appName;

   // qDebug() << msg << endl;

   emit updateProgressSignal(progress, msg.toStdString(), 100);
}
