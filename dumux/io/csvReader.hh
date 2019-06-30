#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
//https://thispointer.com/how-to-read-data-from-a-csv-file-in-c/
/*
 * A class to read data from a csv file.
 */
class CSVReader
{
	std::string delimeter;

public:
	CSVReader(std::string delm = "\t") : delimeter(delm)
	{}
	void getDataFromFile(std::string fileName)
	{
		fileName_=fileName;
		std::ifstream file(fileName_);
		if (file.fail())
			DUNE_THROW(Dune::IOError, "Cannot open file " << fileName_ << "\n");
		int lineNo = 0;
		// std::cout<<"Creat datalist.."<<"\n";
		std::string line = "";
		// Iterate through each line and split the content using delimeter
		while (getline(file, line))
		{
			//std::cout<< "line "<<line<<"\n";
			if (lineNo == 0)
			{
				boost::algorithm::split(headerNames_, line, boost::is_any_of(delimeter));
			}
			else
			{
				std::vector<std::string> vec;
				boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
				std::vector<double> dv;
	            std::transform(vec.begin(), vec.end(), back_inserter(dv), [](const std::string & astr){ return stod( astr) ; } );
				valueArray_.push_back(dv);
			}
			lineNo ++;
		}
		// Close the File
		file.close();
		totalStep_ = lineNo-1;
	}
	// Function to fetch data from a CSV File
	/*
	* Parses through csv file line by line and returns the data
	* in vector of vector of strings.
	*/
	std::vector<std::vector<std::string> > getData()
	{
		std::ifstream file(fileName_);
		if (file.fail())
			DUNE_THROW(Dune::IOError, "Cannot open file " << fileName_ << "\n");
		std::vector<std::vector<std::string> > dataList;
		// std::cout<<"Creat datalist.."<<"\n";
		std::string line = "";
		// Iterate through each line and split the content using delimeter
		while (getline(file, line))
		{
			std::vector<std::string> vec;
			boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
			//std::cout<<"Get line.."<<"\n";
			dataList.push_back(vec);
		}
		// Close the File
		file.close();
		return dataList;
	}

	std::vector<std::vector<double> > getAllValue()
	{
		return valueArray_;
	}

	int getTotalSteps() const
	{
		return totalStep_;
	}
	double getValueAtTime(const std::string &headerName, const double &time) const
	{
		int columeIdx = getIndex(headerName);
		double scanningTime = 0;
		bool notFound = true;
		double lookUpValue = -999;
		int timeStepIdx = 0;
		//std::cout<< time<<" "<<headerName<<" "<<lookUpValue<<" "<<stepHeaderIndex_<<"\n";
		while ((notFound))
		{
			//std::cout<< " "<<scanningTime<<" "<<timeStepIdx<<" "<<totalStep_-1<<"\n";
			scanningTime = valueArray_[timeStepIdx][stepHeaderIndex_];
			if (scanningTime > time)
			{
				notFound = false;
				if (timeStepIdx == 0)
					lookUpValue = valueArray_[timeStepIdx][columeIdx];
				else
					lookUpValue = valueArray_[timeStepIdx-1][columeIdx];
			}
			else
			if (timeStepIdx == totalStep_-1)
			{
				notFound = false;
				lookUpValue = valueArray_[timeStepIdx][columeIdx];
			}
			timeStepIdx++;
		}
		//std::cout<<"DONE "<< time<<" "<<headerName<<" "<<lookUpValue<<"\n";
		return lookUpValue;
	}
	//double getValueAtTime(const std::string &headerName, const double &time) const
	//{
	//	int columeIdx = getIndex(headerName);
	//	double scanningTime = 0;
	//	bool lookUp = true;
	//	double lookUpValue = -999;
	//	int timeStepIdx = 0;
//
	//	while ((lookUp)and(timeStepIdx<totalStep_))
	//	{
	//		scanningTime = valueArray_[timeStepIdx][stepHeaderIndex_];
	//		if (scanningTime <= time) lookUpValue = valueArray_[timeStepIdx][columeIdx];
	//		else lookUp = false;
	//		timeStepIdx++;
	//	}
	//	if (timeStepIdx > totalStep_)
	//		DUNE_THROW(Dune::IOError, "Cannot find value of """ << headerName <<""" at time "<< time << "\n");
//
	//	return lookUpValue;
	//}

	std::vector<double> getValueInHeader(const std::string &headerName) const
	{
		int columeIdx = getIndex(headerName);
		std::vector<double> getValueInHeader;
		for (int stepIdx = 0; stepIdx < totalStep_; stepIdx++)
			getValueInHeader.push_back(valueArray_[stepIdx][columeIdx]);
		return getValueInHeader;
	}

	int getDepthIndex(const double &depth) const
	{
		int depthIdx = 0;
		if (depth <= valueArray_[depthIdx][stepHeaderIndex_]) return depthIdx;
		while ((depthIdx < totalStep_) and (depth >  valueArray_[depthIdx][stepHeaderIndex_]))
			depthIdx++;
		return depthIdx-1;
	}

	double getValueAtDepth(const std::string &headerName, const double &depth) const
	{
		return getValueAtTime(headerName, depth);
	}

	double getMaxValueOverTime(const std::string &headerName) const
	{
		int columeIdx = getIndex(headerName);
		double lookUpValue = -999;
		int timeStepIdx = 0;

		while (timeStepIdx<totalStep_)
		{
			//std::cout<< valueArray_[timeStepIdx][columeIdx] << "\n";
			if (lookUpValue < valueArray_[timeStepIdx][columeIdx])
				lookUpValue = valueArray_[timeStepIdx][columeIdx];
			timeStepIdx++;
		}
		return lookUpValue;
	}

	int getIndex(const std::string &headerName) const
	{
		int index = 0;
		while ((index < headerNames_.size()) and (headerNames_[index]!=headerName))
		{
			//std::cout<< "header: "<<headerNames_[index]<< "\n";
			index++;
		}
		//std::cout<< "headerName: "<<headerName<<" "<<index<< "\n";
		//std::cin.ignore();

		if (index > headerNames_.size())
			DUNE_THROW(Dune::IOError, "Cannot find """ << headerName <<""" in the Header field \n");

		return index;
	}

	void setTimeHeader(std::string timeHeader)
	{
		stepHeaderIndex_ = getIndex(timeHeader);
	}
	void setDepthHeader(std::string depthHeader)
	{
		stepHeaderIndex_ = getIndex(depthHeader);
	}
private:
	std::string fileName_;
	std::vector<std::vector<double> > valueArray_;
	int stepHeaderIndex_;
	std::vector<std::string> headerNames_;
	int totalStep_;
};