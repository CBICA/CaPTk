#include "ThreadedExtraction.h"

void ThreadedExtraction::run() {
   ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("true").toString());
   ApplicationPreferences::GetInstance()->SerializePreferences();
   ApplicationPreferences::GetInstance()->DisplayPreferences();

   if (QFile::exists(this->fullPath))
   {
      QZipReader zr(this->fullPath);
      bool ret = zr.extractAll(this->extractPath);

      if(ret)
      {
         //after extraction remove the zip
         bool successfullyremoved = QFile::remove(this->fullPath.toStdString().c_str());
      }
   }

   qDebug() << "Extraction done in background" << this->fullPath << endl;

   emit resultReady(this->appName);
}	