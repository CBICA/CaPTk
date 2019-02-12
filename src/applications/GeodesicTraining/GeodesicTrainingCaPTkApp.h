#ifndef H_CBICA_GEODESIC_TRAINING_CAPTK_APP
#define H_CBICA_GEODESIC_TRAINING_CAPTK_APP

#include <QString>
#include <thread>
#include <string>

#include "GeodesicTrainingSegmentation.h"

#ifdef APP_BASE_CAPTK_H
#include "ApplicationBase.h"
#endif

template<unsigned int Dimensions = 3>
class GeodesicSegmentationCaPTkApp :
#ifdef APP_BASE_CAPTK_H
    public ApplicationBase, 
#endif
    public GeodesicTrainingSegmentation<float, Dimensions>
{
public:
    typedef itk::Image<PixelType, Dimensions> InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef itk::Image<int, Dimensions>       LabelsImageType;
    typedef typename LabelsImageType::Pointer LabelsImagePointer;

    static_assert((Dimensions == 2 || Dimensions == 3), "2D or 3D Images supported");

    explicit GeodesicSegmentationCaPTkApp() 
    {
        connect(this, &GeodesicSegmentationCaPTkApp<Dimensions>::GeodesicTrainingFinished,
                this, &GeodesicSegmentationCaPTkApp<PixelType, Dimensions>::WorkerFinished
        );
    }

	virtual ~GeodesicSegmentationCaPTkApp() {}

    /** Executes the algorithm in a background thread */
    void Run(std::vector<InputImagePointer> inputImages, LabelsImagePointer labelsImage)
    {
        this->Cancel(); // In case it was already running
        m_Thread = std::thread(&worker, this, inputImages, labelsImage);
    }

    /** Cancels the background thread if it's running */
    void Cancel()
    {
        if (m_Thread.joinable()) {
            m_Thread.join();
        }
    }

public slots:
    /** Just joins the thread (from the main thread). Note: the signal is from this class too. */
    WorkerFinished()
    {
        this->Cancel();
    }

signals:
    GeodesicTrainingFinished(LabelsImagePointer result);
    GeodesicTrainingFinishedWithError(const QString errorMessage);  

private:
    std::thread m_Thread;

    /** The actual background thread */
    void worker(std::vector<InputImagePointer> inputImages, LabelsImagePointer labelsImage)
    {
        this->SetInputImages(inputImages);
        this->SetLabels(labelsImage);
        auto executeResult = this->Execute();

        if (executeResult->ok) {
            emit GeodesicTrainingFinished( executeResult->labelsImage ); // The result segmantation
        }
        else {
            emit GeodesicTrainingFinishedWithError( QString(executeResult->errorMessage) );
        }
    }

    /** Overriden from GeodesicTrainingSegmentation class */
    void progressUpdate(std::string message, int progress) override
    {
#ifdef APP_BASE_CAPTK_H
        messageUpdate(QString::fromStdString(message));
        progressUpdate(progress);
#endif
    }

};

#endif // ! H_CBICA_GEODESIC_TRAINING_CAPTK_APP