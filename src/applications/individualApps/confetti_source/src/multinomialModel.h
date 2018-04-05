/**
\file  multinomialModel.h

\brief Implementation of the Confetti Algorithm

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2017 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html
Author: Ratheesh Kalarot
*/
#include "dtiUtils.h"
using namespace cv;
/**
\class Confetti
\brief Implementation of  Multinomial Model part of Confetti algorithm
Reference:
Tunç B, Smith AR, Wasserman D, et al.
Multinomial Probabilistic Fiber Representation for Connectivity Driven Clustering.
Information processing in medical imaging?: proceedings of the . conference.
2013;23:730-741.
*/
class MNM
{
public:
  //Set listner for progress update
  void setListner(TXT_CALLBACK_TYPE listner)
  {
    Utils::setListner(listner);
    return;
  }
  //gen Adaptive Cluster
  int genAdaptiveCluster(string signatureFile, string tempRoot, string outFile, int K = 200, bool bUseTemplate = false, string outTemplateDir = std::string())
	{
		Utils::resetTimer();
		Utils::logText("Reading:" + signatureFile);
		vector< vector<double> > data = Utils::readCsv<double>(signatureFile);
		if (data.size() == 0)
		{
			Utils::logText("Error reading:" + signatureFile);
			return EXIT_FAILURE;
		}
		// remove fibers with a few possible connections
		Utils::logText("Pruning fibers with a few possible connections");
		vector< vector<double> > dataPruned;
		vector<int> selectedIds;
		for (size_t r = 0; r < data.size(); r++)
		{
			if (*max_element(data[r].begin(), data[r].end()) > 50)
			{
				dataPruned.push_back(data[r]);
				selectedIds.push_back(r);
			}
		}
		size_t realNumFib = data.size();
		Utils::logText(Utils::toString(data.size() - dataPruned.size()) + " fibers were removed out of" + Utils::toString(data.size()));
		data = dataPruned;
		Mat dataMat = Utils::toMat(data);

		// transform data to facilitate robust clustering
		Utils::logText("Applying DfIdf Transform");
		Utils::dfIdfTransform(dataMat);
		Mat mx;
		cv::reduce(dataMat, mx, 1, CV_REDUCE_MAX);
		for (int r = 0; r < mx.rows; r++)
		{
			dataMat.row(r) = dataMat.row(r) * 100 / mx.at<double>(r,0);
		}
		transform(dataMat.begin<double>(), dataMat.end<double>(), dataMat.begin<double>(), floorf);
    Mat prior_alphas, prior_priors;
    if (bUseTemplate)
    {
      prior_alphas = Utils::toMat(Utils::readCsv<double>(tempRoot + "/template_Alp.csv"));
      prior_priors = Utils::toMat(Utils::readCsv<double>(tempRoot + "/template_Pr.csv"));
      K = prior_priors.rows;
    }


		Utils::logText("Running Optimization");
    vector<int> clusterID_;
    if (bUseTemplate)
    {
      clusterID_ = fitMixture(dataMat, K, 500, 1e-3, prior_alphas, prior_priors, true);
    }
    else
    {
      clusterID_ = fitMixture(dataMat, K, 500, 1e-3, prior_alphas, prior_priors, false);
    }
   
		vector<int> clusterID = vector<int>(realNumFib, -1);
		for (size_t idx = 0; idx < clusterID_.size(); idx++)
		{
			clusterID[selectedIds[idx]] = clusterID_[idx];
		}
		Utils::logText("Writing output:" + outFile);
		Utils::writeCsv<int>(clusterID, outFile);
		Utils::logText("Finished adaptive clustering");
    if (!outTemplateDir.empty())
    {
      auto alphas = Utils::toVector<double>(Mixture_alpha);
      auto priors = Utils::toVector<double>(Mixture_prior);
      Utils::writeCsv(alphas, outTemplateDir + "/template_Alp.csv");
      Utils::writeCsv(priors, outTemplateDir + "/template_Pr.csv");
    }

    return EXIT_SUCCESS;
	}

private:

  Mat logsum(const Mat &dataIn)
  {
    //dataIn ->200*512 
    Mat dataTemp = dataIn.clone();
    Mat amax; //-> 512
    cv::reduce(dataTemp, amax, 0, CV_REDUCE_MAX);
    for (int r = 0; r < dataTemp.rows; r++)
    {
      dataTemp.row(r) = dataTemp.row(r) - amax;
    }
    Mat datExp, datExpSum, datLog;
    cv::exp(dataTemp, datExp);
    cv::reduce(datExp, datExpSum, 0, CV_REDUCE_SUM);
    cv::log(datExpSum, datLog);
    Mat s = amax + datLog;
    //double sum3 = sum(s).val[0];
    return s;
  }
  Mat evaluate(const Mat &Data, double& sumOut)
  {
    Mat n;
    cv::reduce(Data, n, 1, CV_REDUCE_SUM);
    Mat logMat;
    cv::log(Mixture_alpha.t(), logMat);
    Mat term1 = Data* logMat;
    Mat sumMat;
    cv::reduce(Utils::gammaLn(Data + 1), sumMat, 1, CV_REDUCE_SUM);
    Mat part1 = Utils::gammaLn(n + 1);
    Mat term2 = part1 - sumMat;
    Mat logL = term1.clone();
    for (int c = 0; c < logL.cols; c++)
    {
      logL.col(c) = logL.col(c) + term2;
    }
    logL = logL.t();// rows=200, cols=512
    Mat logP_;
    Utils::threshTruncVal(Mixture_prior, 0, DBL_EPSILON);
    cv::log(Mixture_prior, logP_);//rows=200, cols=1
    Mat logP = logL.clone();
    for (int c = 0; c < logL.cols; c++)//logP = logL + logP rows=200 cols 512
    {
      logP.col(c) = logP.col(c) + logP_;
    }
    Mat logLike = logsum(logP);//->512
    Mat R;//->512,200
    Mat logPNew = logP.clone();
    for (int r = 0; r < logP.rows; r++)
    {
      logPNew.row(r) = logP.row(r) - logLike;
    }
    cv::exp(logPNew, R);

    sumOut = cv::sum(logLike).val[0];
    Mat NkDbg;
    cv::reduce(R, NkDbg, 0, CV_REDUCE_SUM);
    NkDbg = NkDbg + 1e-12;
    return R.t();
  }
  vector<int> fitMixture(Mat Data, int K = 200, int maxIt = 100, double threshold = 1e-3, Mat Prior_alpha = Mat(), Mat Prior_priors = Mat(), bool adaptive = false)
  {
    //int K;
    Mat A;
    int numD = Data.rows;
    int numF = Data.cols;
    if (!Prior_alpha.empty())
    {
      K = Prior_alpha.cols;
      Mixture_prior = Prior_priors.clone();
      Mixture_alpha = Prior_alpha.clone() + 2 * 1e-12;
      A = Prior_alpha;
      Utils::logText("Using prior knowledge");
    }
    else
    {
      Mat  temp(K, numF, CV_64FC1); // Or: Mat mat(2, 4, CV_64FC1);
      randu(temp, Scalar(0.0), Scalar(1.0));
      Mixture_alpha = temp;
      Mixture_prior = Mat(vector<double>(K, 1.0 / K));
    }
    double convergence = 7.0;
    int it = 0;
    double logLike0;
    Mat R = evaluate(Data, logLike0);
    Mat dsum;
    cv::reduce(Data, dsum, 1, CV_REDUCE_SUM);
    while (it<maxIt && convergence>threshold)
    {
      Utils::progressUpdate(it + 1, maxIt, "Run Optimization");
      //E step calculate priors
      Mat Nk;//->200L
      cv::reduce(R, Nk, 0, CV_REDUCE_SUM);
      Nk = Nk + 1e-12;
      double sumNK = cv::sum(Nk).val[0];
      Nk = Nk * (numD / sumNK);

      Mat prior = Nk.t() / double(numD);
      Utils::threshTruncVal(prior, 0, 1e-12);//prior[prior <= 0] = 1e-12
      prior = prior / sum(prior).val[0];


      //M step calculate alphas
      Mat term1 = R.t()* Data; //200Lx87L
      Mat term2 = R.t()* dsum + 1e-12; //200L
      if (adaptive)
      {
        Mat pterm1 = 5000000 * Prior_alpha;
        Mat pterm2;
        cv::reduce(pterm1, pterm2, 1, CV_REDUCE_SUM);
        term1 = term1 + pterm1;
        term2 = term2 + pterm2;
      }
      Mat newA = term1;
      for (int c = 0; c < newA.cols; c++)
      {
        newA.col(c) = newA.col(c) / term2;
      }
      Mixture_prior = prior;
      Mixture_alpha = newA + 1e-12;

      //Convergence
      double logLike1;
      R = evaluate(Data, logLike1);
      convergence = abs(logLike1 - logLike0);
      logLike0 = logLike1;
      it += 1;
      Utils::logText("it:" + Utils::toString(it) + " logL:" + Utils::toString(logLike0) + " Diff:" + Utils::toString(convergence));
    }
    //R 512 * 200 lables= 512
    vector<int> labels;
    for (int r = 0; r < R.rows; r++)
    {
      double a, b;
      cv::Point minLoc, maxLoc;
      cv::minMaxLoc(R.row(r), &a, &b, &minLoc, &maxLoc);
      labels.push_back(max(maxLoc.x, maxLoc.y));
    }
    // calculate final mixture using only data
    if (adaptive)
    {
      Mat term1 = R.t()* Data; //200Lx87L
      Mat term2 = R.t()* dsum + 1e-12; //200L
      Mat newA = term1;
      for (int c = 0; c < newA.cols; c++)
      {
        newA.col(c) = newA.col(c) / term2;
      }
      //Mixture_prior = prior;
      Mixture_alpha = newA + 1e-12;
    }
    return labels;
  }

  //member variables
	Mat Mixture_prior;
	Mat Mixture_alpha;

};
