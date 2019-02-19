#include "GeodesicTrainingCaPTkApp.h"

GeodesicTrainingApplicationBase::GeodesicTrainingApplicationBase(QObject* parent) : QObject(parent)
{
	connect(this, SIGNAL(GeodesicTrainingFinished()),
		this, SLOT(WorkerFinished())
	);

	connect(this, SIGNAL(GeodesicTrainingFinishedWithError()),
		this, SLOT(WorkerFinished())
	);
}

void GeodesicTrainingApplicationBase::Cancel()
{
	if (m_Thread.joinable()) {
		m_Thread.join();
	}
}

void GeodesicTrainingApplicationBase::WorkerFinished()
{
	this->Cancel();
}

// GeodesicSegmentationCaPTkApp is a templated class