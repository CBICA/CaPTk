#ifndef H_CBICA_GEODESIC_TRAINING_CAPTK_APP
#define H_CBICA_GEODESIC_TRAINING_CAPTK_APP

#include <QObject>
#include <QString>
#include <thread>
#include <string>

#include "GeodesicTrainingSegmentation.h"

//#include "ApplicationBase.h"
#include <QObject>

class GeodesicTrainingApplicationBase : public QObject
{
	Q_OBJECT

public:
	GeodesicTrainingApplicationBase(QObject* parent = nullptr);

	/** Cancels the background thread if it's running */
	virtual void Cancel();

public slots:
	/** Just joins the thread (from the main thread). Note: the signal is from this class too. */
	void WorkerFinished();

signals:
	void GeodesicTrainingFinished();
	//void GeodesicTrainingFinished3D(typename itk::Image<int, 3>::Pointer result);
	//void GeodesicTrainingFinished2D(typename itk::Image<int, 2>::Pointer result);
	void GeodesicTrainingFinishedWithError(QString errorMessage);
	void GeodesicTrainingProgressUpdate(int progress, std::string message, int max);

protected:
	std::thread m_Thread;
};

template<unsigned int Dimensions = 3>
class GeodesicTrainingCaPTkApp :
	public GeodesicTrainingApplicationBase,
	public GeodesicTrainingSegmentation::Coordinator<float, Dimensions>
{
public:
    typedef itk::Image<float, Dimensions>     InputImageType;
    typedef typename InputImageType::Pointer  InputImagePointer;
    typedef itk::Image<int, Dimensions>       LabelsImageType;
    typedef typename LabelsImageType::Pointer LabelsImagePointer;

    explicit GeodesicTrainingCaPTkApp(QObject* parent = nullptr) : GeodesicTrainingApplicationBase(parent) {}

	virtual ~GeodesicTrainingCaPTkApp() {}

    /** Executes the algorithm in a background thread */
    void Run(std::vector<InputImagePointer> inputImages, LabelsImagePointer labelsImage)
    {
        this->Cancel(); // In case it was already running
		m_Thread = std::thread(&GeodesicTrainingCaPTkApp<Dimensions>::worker, this, inputImages, labelsImage);
    }

	/** Overriden from GeodesicTrainingSegmentation class */
	void progressUpdate(std::string message, int progress) override
	{
		emit GeodesicTrainingProgressUpdate(progress, message, 100);
	}

private:

    /** The actual background thread */
    void worker(std::vector<InputImagePointer> inputImages, LabelsImagePointer labelsImage)
    {
        this->SetInputImages(inputImages);
        this->SetLabels(labelsImage);
		this->SetSaveAll(true);
		this->SetProcessing(true, 6, true, false, false, 5000000); // Disable pixel limit, all else default.
        auto executeResult = this->Execute();

        if (executeResult->ok) {
            /*if (Dimensions == 3) { emit GeodesicTrainingFinished3D( executeResult->labelsImage ); }
            else {                 emit GeodesicTrainingFinished2D( executeResult->labelsImage ); }*/
        
			emit GeodesicTrainingProgressUpdate(100, "Geodesic Training Segmentation: Finished!", 100);
			emit GeodesicTrainingFinished();
		}
        else {
            emit GeodesicTrainingFinishedWithError( QString(executeResult->errorMessage.c_str()) );			
        }
    }

};

#endif // ! H_CBICA_GEODESIC_TRAINING_CAPTK_APP