/*
Copyright (c) 2010-2016, Mathieu Labbe - IntRoLab - Universite de Sherbrooke
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Universite de Sherbrooke nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <rtabmap/core/FlannIndex.h>
#include <rtabmap/utilite/ULogger.h>

#include "rtflann/flann.hpp"
#include <sys/mman.h>

namespace rtabmap {

FlannIndex::FlannIndex():
		index_(0),
		nextIndex_(0),
		featuresType_(0),
		featuresDim_(0),
		isLSH_(false),
		useDistanceL1_(false),
		rebalancingFactor_(2.0f)
{
}
FlannIndex::~FlannIndex()
{
	this->release();
}

void FlannIndex::release()
{
	UDEBUG("");
	if(index_)
	{
		if(featuresType_ == CV_8UC1)
		{
			delete (rtflann::Index<rtflann::Hamming<unsigned char> >*)index_;
		}
		else
		{
			if(useDistanceL1_)
			{
				delete (rtflann::Index<rtflann::L1<float> >*)index_;
			}
			else if(featuresDim_ <= 3)
			{
				delete (rtflann::Index<rtflann::L2_Simple<float> >*)index_;
			}
			else
			{
				delete (rtflann::Index<rtflann::L2<float> >*)index_;
			}
		}
		index_ = 0;
	}
	nextIndex_ = 0;
	isLSH_ = false;
	addedDescriptors_.clear();
	removedIndexes_.clear();
	UDEBUG("");
}

void FlannIndex::debug()
{
	std::cout << "flann size: " << ((const rtflann::Index<rtflann::L1<float> >*)index_)->size() << std::endl;
	((rtflann::Index<rtflann::L1<float> >*)index_)->debug();
	std::cout << "size of addedDescriptors_: " << addedDescriptors_.size() << std::endl;
	std::cout << "size of removedIndexes_: " << removedIndexes_.size() << std::endl;
}

void FlannIndex::save_index()
{
	std::ofstream *outfile;
	outfile = new std::ofstream();
	outfile->open("flann.dat", std::ios::out | std::ios::binary | std::ios::trunc);
	((rtflann::Index<rtflann::L1<float>> *)index_)->save_index(outfile);

	//FIXME - save variables in FlannIndex
	//unsigned int nextIndex_;
	outfile->write(reinterpret_cast<char *>(&nextIndex_), sizeof(unsigned int));

	//int featuresType_;
	outfile->write(reinterpret_cast<char *>(&featuresType_), sizeof(int));

	// int featuresDim_;
	outfile->write(reinterpret_cast<char *>(&featuresDim_), sizeof(int));

	// bool isLSH_;
	outfile->write(reinterpret_cast<char *>(&isLSH_), sizeof(bool));

	// bool useDistanceL1_; 
	outfile->write(reinterpret_cast<char *>(&useDistanceL1_), sizeof(bool));

	// float rebalancingFactor_;
	outfile->write(reinterpret_cast<char *>(&rebalancingFactor_), sizeof(float));

	// // keep feature in memory until the tree is rebuilt
	// // (in case the word is deleted when removed from the VWDictionary)
	// std::map<int, cv::Mat> addedDescriptors_;
	std::cout << "size of addedDescriptors_" << addedDescriptors_.size() << std::endl;

	// std::list<int> removedIndexes_;

	outfile->close();
}

void FlannIndex::load_index()
{
	std::ifstream *infile;
	infile = new std::ifstream();
	infile->open("/home/jseng/herbie_rtabmap/build/flann.dat", std::ios::in | std::ios::binary);

	//Call mmap() to have the memory block allocated at the same starting address.
	char *addr;
	unsigned long long int starting_addr = 0x600000000000;
	unsigned int length = 1300000000;
	addr = (char *)mmap((void *)starting_addr, length, PROT_READ | PROT_WRITE, MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (addr == MAP_FAILED)
	{
		std::cout << "Error" << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		std::cout << "Successful mapping to: " << (unsigned long long)addr << std::endl;
	}

	//read in all the visual words from a file
	int point_size;
	infile->read(reinterpret_cast<char *>(&point_size), sizeof(int));

	float *data_ptr = reinterpret_cast<float *>(addr);
	for (int i = 0; i < point_size; i++)
	{
		infile->read(reinterpret_cast<char *>(data_ptr), 256 * sizeof(float));
		data_ptr += 256;
	}

	//create a Matrix from the data
	cv::Mat features = cv::Mat(point_size, 256, CV_32F, addr);
	featuresType_ = features.type();
	featuresDim_ = features.cols;
	useDistanceL1_ = true;
	rebalancingFactor_ = 2.0;

	rtflann::KDTreeIndexParams params(4);

	rtflann::Matrix<float> dataset((float*)features.data, features.rows, features.cols);
	index_ = new rtflann::Index<rtflann::L1<float> >(dataset, params);

	((rtflann::Index<rtflann::L1<float>> *)index_)->set_cached(1); //indicate to use a cache index
	((rtflann::Index<rtflann::L1<float>> *)index_)->buildIndex();

	char* var_start_ptr = (char*) (starting_addr + point_size*256*sizeof(float)); 
	((rtflann::Index<rtflann::L1<float>> *)index_)->load_index(infile, var_start_ptr);

	//FIXME - load variables in FlannIndex
	//unsigned int nextIndex_;
	infile->read(reinterpret_cast<char *>(&nextIndex_), sizeof(unsigned int));

	//int featuresType_;
	infile->read(reinterpret_cast<char *>(&featuresType_), sizeof(int));

	// int featuresDim_;
	infile->read(reinterpret_cast<char *>(&featuresDim_), sizeof(int));

	// bool isLSH_;
	infile->read(reinterpret_cast<char *>(&isLSH_), sizeof(bool));

	// bool useDistanceL1_;
	infile->read(reinterpret_cast<char *>(&useDistanceL1_), sizeof(bool));

	// float rebalancingFactor_;
	infile->read(reinterpret_cast<char *>(&rebalancingFactor_), sizeof(float));

	infile->close();
}

unsigned int FlannIndex::indexedFeatures() const
{
	if(!index_)
	{
		return 0;
	}
	if(featuresType_ == CV_8UC1)
	{
		return ((const rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->size();
	}
	else
	{
		if(useDistanceL1_)
		{
			return ((const rtflann::Index<rtflann::L1<float> >*)index_)->size();
		}
		else if(featuresDim_ <= 3)
		{
			return ((const rtflann::Index<rtflann::L2_Simple<float> >*)index_)->size();
		}
		else
		{
			return ((const rtflann::Index<rtflann::L2<float> >*)index_)->size();
		}
	}
}

// return Bytes
unsigned long FlannIndex::memoryUsed() const
{
	if(!index_)
	{
		return 0;
	}
	unsigned long memoryUsage = sizeof(FlannIndex);
	memoryUsage += addedDescriptors_.size() * (sizeof(int) + sizeof(cv::Mat) + sizeof(std::map<int, cv::Mat>::iterator)) + sizeof(std::map<int, cv::Mat>);
	memoryUsage += sizeof(std::list<int>) + removedIndexes_.size() * sizeof(int);
	if(featuresType_ == CV_8UC1)
	{
		memoryUsage += ((const rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->usedMemory();
	}
	else
	{
		if(useDistanceL1_)
		{
			memoryUsage += ((const rtflann::Index<rtflann::L1<float> >*)index_)->usedMemory();
		}
		else if(featuresDim_ <= 3)
		{
			memoryUsage += ((const rtflann::Index<rtflann::L2_Simple<float> >*)index_)->usedMemory();
		}
		else
		{
			memoryUsage += ((const rtflann::Index<rtflann::L2<float> >*)index_)->usedMemory();
		}
	}
	return memoryUsage;
}

void FlannIndex::buildLinearIndex(
		const cv::Mat & features,
		bool useDistanceL1,
		float rebalancingFactor)
{
	UDEBUG("");
	this->release();
	UASSERT(index_ == 0);
	UASSERT(features.type() == CV_32FC1 || features.type() == CV_8UC1);
	featuresType_ = features.type();
	featuresDim_ = features.cols;
	useDistanceL1_ = useDistanceL1;
	rebalancingFactor_ = rebalancingFactor;

	rtflann::LinearIndexParams params;

	if(featuresType_ == CV_8UC1)
	{
		rtflann::Matrix<unsigned char> dataset(features.data, features.rows, features.cols);
		index_ = new rtflann::Index<rtflann::Hamming<unsigned char> >(dataset, params);
		((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->buildIndex();
	}
	else
	{
		rtflann::Matrix<float> dataset((float*)features.data, features.rows, features.cols);
		if(useDistanceL1_)
		{
			index_ = new rtflann::Index<rtflann::L1<float> >(dataset, params);
			((rtflann::Index<rtflann::L1<float> >*)index_)->buildIndex();
		}
		else if(featuresDim_ <=3)
		{
			index_ = new rtflann::Index<rtflann::L2_Simple<float> >(dataset, params);
			((rtflann::Index<rtflann::L2_Simple<float> >*)index_)->buildIndex();
		}
		else
		{
			index_ = new rtflann::Index<rtflann::L2<float> >(dataset, params);
			((rtflann::Index<rtflann::L2<float> >*)index_)->buildIndex();
		}
	}

	// incremental FLANN: we should add all headers separately in case we remove
	// some indexes (to keep underlying matrix data allocated)
	if(rebalancingFactor_ > 1.0f)
	{
		for(int i=0; i<features.rows; ++i)
		{
			addedDescriptors_.insert(std::make_pair(nextIndex_++, features.row(i)));
		}
	}
	else
	{
		// tree won't ever be rebalanced, so just keep only one header for the data
		addedDescriptors_.insert(std::make_pair(nextIndex_, features));
		nextIndex_ += features.rows;
	}
	UDEBUG("");
}

void FlannIndex::buildKDTreeIndex(
		const cv::Mat & features,
		int trees,
		bool useDistanceL1,
		float rebalancingFactor)
{
	UDEBUG("");
	this->release();
	UASSERT(index_ == 0);
	UASSERT(features.type() == CV_32FC1 || features.type() == CV_8UC1);
	featuresType_ = features.type();
	featuresDim_ = features.cols;
	useDistanceL1_ = useDistanceL1;
	rebalancingFactor_ = rebalancingFactor;

	rtflann::KDTreeIndexParams params(trees);

	if(featuresType_ == CV_8UC1)
	{
		rtflann::Matrix<unsigned char> dataset(features.data, features.rows, features.cols);
		index_ = new rtflann::Index<rtflann::Hamming<unsigned char> >(dataset, params);
		((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->buildIndex();
	}
	else
	{
		rtflann::Matrix<float> dataset((float*)features.data, features.rows, features.cols);
		if(useDistanceL1_)
		{
			index_ = new rtflann::Index<rtflann::L1<float> >(dataset, params);
			((rtflann::Index<rtflann::L1<float> >*)index_)->buildIndex();
		}
		else if(featuresDim_ <=3)
		{
			index_ = new rtflann::Index<rtflann::L2_Simple<float> >(dataset, params);
			((rtflann::Index<rtflann::L2_Simple<float> >*)index_)->buildIndex();
		}
		else
		{
			index_ = new rtflann::Index<rtflann::L2<float> >(dataset, params);
			((rtflann::Index<rtflann::L2<float> >*)index_)->buildIndex();
		}
	}

	// incremental FLANN: we should add all headers separately in case we remove
	// some indexes (to keep underlying matrix data allocated)
	if(rebalancingFactor_ > 1.0f)
	{
		for(int i=0; i<features.rows; ++i)
		{
			addedDescriptors_.insert(std::make_pair(nextIndex_++, features.row(i)));
		}
	}
	else
	{
		// tree won't ever be rebalanced, so just keep only one header for the data
		addedDescriptors_.insert(std::make_pair(nextIndex_, features));
		nextIndex_ += features.rows;
	}
	UDEBUG("");
}

void FlannIndex::buildKDTreeSingleIndex(
		const cv::Mat & features,
		int leafMaxSize,
		bool reorder,
		bool useDistanceL1,
		float rebalancingFactor)
{
	UDEBUG("");
	this->release();
	UASSERT(index_ == 0);
	UASSERT(features.type() == CV_32FC1 || features.type() == CV_8UC1);
	featuresType_ = features.type();
	featuresDim_ = features.cols;
	useDistanceL1_ = useDistanceL1;
	rebalancingFactor_ = rebalancingFactor;

	rtflann::KDTreeSingleIndexParams params(leafMaxSize, reorder);

	if(featuresType_ == CV_8UC1)
	{
		rtflann::Matrix<unsigned char> dataset(features.data, features.rows, features.cols);
		index_ = new rtflann::Index<rtflann::Hamming<unsigned char> >(dataset, params);
		((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->buildIndex();
	}
	else
	{
		rtflann::Matrix<float> dataset((float*)features.data, features.rows, features.cols);
		if(useDistanceL1_)
		{
			index_ = new rtflann::Index<rtflann::L1<float> >(dataset, params);
			((rtflann::Index<rtflann::L1<float> >*)index_)->buildIndex();
		}
		else if(featuresDim_ <=3)
		{
			index_ = new rtflann::Index<rtflann::L2_Simple<float> >(dataset, params);
			((rtflann::Index<rtflann::L2_Simple<float> >*)index_)->buildIndex();
		}
		else
		{
			index_ = new rtflann::Index<rtflann::L2<float> >(dataset, params);
			((rtflann::Index<rtflann::L2<float> >*)index_)->buildIndex();
		}
	}

	// incremental FLANN: we should add all headers separately in case we remove
	// some indexes (to keep underlying matrix data allocated)
	if(rebalancingFactor_ > 1.0f)
	{
		for(int i=0; i<features.rows; ++i)
		{
			addedDescriptors_.insert(std::make_pair(nextIndex_++, features.row(i)));
		}
	}
	else
	{
		// tree won't ever be rebalanced, so just keep only one header for the data
		addedDescriptors_.insert(std::make_pair(nextIndex_, features));
		nextIndex_ += features.rows;
	}
	UDEBUG("");
}

void FlannIndex::buildLSHIndex(
		const cv::Mat & features,
		unsigned int table_number,
		unsigned int key_size,
		unsigned int multi_probe_level,
		float rebalancingFactor)
{
	UDEBUG("");
	this->release();
	UASSERT(index_ == 0);
	UASSERT(features.type() == CV_8UC1);
	featuresType_ = features.type();
	featuresDim_ = features.cols;
	useDistanceL1_ = true;
	rebalancingFactor_ = rebalancingFactor;

	rtflann::Matrix<unsigned char> dataset(features.data, features.rows, features.cols);
	index_ = new rtflann::Index<rtflann::Hamming<unsigned char> >(dataset, rtflann::LshIndexParams(12, 20, 2));
	((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->buildIndex();

	// incremental FLANN: we should add all headers separately in case we remove
	// some indexes (to keep underlying matrix data allocated)
	if(rebalancingFactor_ > 1.0f)
	{
		for(int i=0; i<features.rows; ++i)
		{
			addedDescriptors_.insert(std::make_pair(nextIndex_++, features.row(i)));
		}
	}
	else
	{
		// tree won't ever be rebalanced, so just keep only one header for the data
		addedDescriptors_.insert(std::make_pair(nextIndex_, features));
		nextIndex_ += features.rows;
	}
	UDEBUG("");
}

bool FlannIndex::isBuilt()
{
	return index_!=0;
}

std::vector<unsigned int> FlannIndex::addPoints(const cv::Mat & features)
{
	if(!index_)
	{
		UERROR("Flann index not yet created!");
		return std::vector<unsigned int>();
	}
	UASSERT(features.type() == featuresType_);
	UASSERT(features.cols == featuresDim_);
	bool indexRebuilt = false;
	size_t removedPts = 0;
	if(featuresType_ == CV_8UC1)
	{
		rtflann::Matrix<unsigned char> points(features.data, features.rows, features.cols);
		rtflann::Index<rtflann::Hamming<unsigned char> > * index = (rtflann::Index<rtflann::Hamming<unsigned char> >*)index_;
		removedPts = index->removedCount();
		index->addPoints(points, 0);
		// Rebuild index if it is now X times in size
		if(rebalancingFactor_ > 1.0f && size_t(float(index->sizeAtBuild()) * rebalancingFactor_) < index->size()+index->removedCount())
		{
			UDEBUG("Rebuilding FLANN index: %d -> %d", (int)index->sizeAtBuild(), (int)(index->size()+index->removedCount()));
			index->buildIndex();
            std::cout << "---------------byte - addPoints() - FlannIndex.cpp--------------" << std::endl;
		}
		// if no more removed points, the index has been rebuilt
		indexRebuilt = index->removedCount() == 0 && removedPts>0;
	}
	else
	{
		rtflann::Matrix<float> points((float*)features.data, features.rows, features.cols);
		if(useDistanceL1_)
		{
			rtflann::Index<rtflann::L1<float> > * index = (rtflann::Index<rtflann::L1<float> >*)index_;
			removedPts = index->removedCount();
			index->addPoints(points, 0);
			// Rebuild index if it doubles in size
			if(rebalancingFactor_ > 1.0f && size_t(float(index->sizeAtBuild()) * rebalancingFactor_) < index->size()+index->removedCount())
			{
                std::cout << "---------------L1 - addPoints() - FlannIndex.cpp--------------" << std::endl;
				UDEBUG("Rebuilding FLANN index: %d -> %d", (int)index->sizeAtBuild(), (int)(index->size()+index->removedCount()));
				index->buildIndex();
			}
			// if no more removed points, the index has been rebuilt
			indexRebuilt = index->removedCount() == 0 && removedPts>0;
		}
		else if(featuresDim_ <= 3)
		{
			rtflann::Index<rtflann::L2_Simple<float> > * index = (rtflann::Index<rtflann::L2_Simple<float> >*)index_;
			removedPts = index->removedCount();
			index->addPoints(points, 0);
			// Rebuild index if it doubles in size
			if(rebalancingFactor_ > 1.0f && size_t(float(index->sizeAtBuild()) * rebalancingFactor_) < index->size()+index->removedCount())
			{
				UDEBUG("Rebuilding FLANN index: %d -> %d", (int)index->sizeAtBuild(), (int)(index->size()+index->removedCount()));
				index->buildIndex();
                std::cout << "---------------L2_simple - addPoints() - FlannIndex.cpp--------------" << std::endl;
			}
			// if no more removed points, the index has been rebuilt
			indexRebuilt = index->removedCount() == 0 && removedPts>0;
		}
		else
		{
			rtflann::Index<rtflann::L2<float> > * index = (rtflann::Index<rtflann::L2<float> >*)index_;
			removedPts = index->removedCount();
			index->addPoints(points, 0);
			// Rebuild index if it doubles in size
			if(rebalancingFactor_ > 1.0f && size_t(float(index->sizeAtBuild()) * rebalancingFactor_) < index->size()+index->removedCount())
			{
				UDEBUG("Rebuilding FLANN index: %d -> %d", (int)index->sizeAtBuild(), (int)(index->size()+index->removedCount()));
				index->buildIndex();
                std::cout << "---------------L2 - addPoints() - FlannIndex.cpp--------------" << std::endl;
			}
			// if no more removed points, the index has been rebuilt
			indexRebuilt = index->removedCount() == 0 && removedPts>0;
		}
	}

	if(indexRebuilt)
	{
		UASSERT(removedPts == removedIndexes_.size());
		// clean not used features
		for(std::list<int>::iterator iter=removedIndexes_.begin(); iter!=removedIndexes_.end(); ++iter)
		{
			addedDescriptors_.erase(*iter);
		}
		removedIndexes_.clear();
	}

	// incremental FLANN: we should add all headers separately in case we remove
	// some indexes (to keep underlying matrix data allocated)
	std::vector<unsigned int> indexes;
	for(int i=0; i<features.rows; ++i)
	{
		indexes.push_back(nextIndex_);
		addedDescriptors_.insert(std::make_pair(nextIndex_++, features.row(i)));
	}

	return indexes;
}

void FlannIndex::removePoint(unsigned int index)
{
	if(!index_)
	{
		UERROR("Flann index not yet created!");
		return;
	}

	// If a Segmentation fault occurs in removePoint(), verify that you have this fix in your installed "flann/algorithms/nn_index.h":
	// 707 - if (ids_[id]==id) {
	// 707 + if (id < ids_.size() && ids_[id]==id) {
	// ref: https://github.com/mariusmuja/flann/commit/23051820b2314f07cf40ba633a4067782a982ff3#diff-33762b7383f957c2df17301639af5151

	if(featuresType_ == CV_8UC1)
	{
		((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->removePoint(index);
	}
	else if(useDistanceL1_)
	{
		((rtflann::Index<rtflann::L1<float> >*)index_)->removePoint(index);
	}
	else if(featuresDim_ <= 3)
	{
		((rtflann::Index<rtflann::L2_Simple<float> >*)index_)->removePoint(index);
	}
	else
	{
		((rtflann::Index<rtflann::L2<float> >*)index_)->removePoint(index);
	}

	removedIndexes_.push_back(index);
}

void FlannIndex::knnSearch(
		const cv::Mat & query,
		cv::Mat & indices,
		cv::Mat & dists,
		int knn,
		int checks,
		float eps,
		bool sorted) const
{
	if(!index_)
	{
		UERROR("Flann index not yet created!");
		return;
	}
	indices.create(query.rows, knn, sizeof(size_t)==8?CV_64F:CV_32S);
	dists.create(query.rows, knn, featuresType_ == CV_8UC1?CV_32S:CV_32F);

	rtflann::Matrix<size_t> indicesF((size_t*)indices.data, indices.rows, indices.cols);

	rtflann::SearchParams params = rtflann::SearchParams(checks, eps, sorted);

	if(featuresType_ == CV_8UC1)
	{
		rtflann::Matrix<unsigned int> distsF((unsigned int*)dists.data, dists.rows, dists.cols);
		rtflann::Matrix<unsigned char> queryF(query.data, query.rows, query.cols);
		((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->knnSearch(queryF, indicesF, distsF, knn, params);
	}
	else
	{
		rtflann::Matrix<float> distsF((float*)dists.data, dists.rows, dists.cols);
		rtflann::Matrix<float> queryF((float*)query.data, query.rows, query.cols);
		if(useDistanceL1_)
		{
			((rtflann::Index<rtflann::L1<float> >*)index_)->knnSearch(queryF, indicesF, distsF, knn, params);
		}
		else if(featuresDim_ <= 3)
		{
			((rtflann::Index<rtflann::L2_Simple<float> >*)index_)->knnSearch(queryF, indicesF, distsF, knn, params);
		}
		else
		{
			((rtflann::Index<rtflann::L2<float> >*)index_)->knnSearch(queryF, indicesF, distsF, knn, params);
		}
	}
}

void FlannIndex::radiusSearch(
		const cv::Mat & query,
		std::vector<std::vector<size_t> > & indices,
		std::vector<std::vector<float> > & dists,
		float radius,
		int maxNeighbors,
		int checks,
		float eps,
		bool sorted) const
{
	if(!index_)
	{
		UERROR("Flann index not yet created!");
		return;
	}

	rtflann::SearchParams params = rtflann::SearchParams(checks, eps, sorted);
	params.max_neighbors = maxNeighbors<=0?-1:maxNeighbors; // -1 is all in radius

	if(featuresType_ == CV_8UC1)
	{
		std::vector<std::vector<unsigned int> > distsF;
		rtflann::Matrix<unsigned char> queryF(query.data, query.rows, query.cols);
		((rtflann::Index<rtflann::Hamming<unsigned char> >*)index_)->radiusSearch(queryF, indices, distsF, radius*radius, params);
		dists.resize(distsF.size());
		for(unsigned int i=0; i<dists.size(); ++i)
		{
			dists[i].resize(distsF[i].size());
			for(unsigned int j=0; j<distsF[i].size(); ++j)
			{
				dists[i][j] = (float)distsF[i][j];
			}
		}
	}
	else
	{
		rtflann::Matrix<float> queryF((float*)query.data, query.rows, query.cols);
		if(useDistanceL1_)
		{
			((rtflann::Index<rtflann::L1<float> >*)index_)->radiusSearch(queryF, indices, dists, radius*radius, params);
		}
		else if(featuresDim_ <= 3)
		{
			((rtflann::Index<rtflann::L2_Simple<float> >*)index_)->radiusSearch(queryF, indices, dists, radius*radius, params);
		}
		else
		{
			((rtflann::Index<rtflann::L2<float> >*)index_)->radiusSearch(queryF, indices, dists, radius*radius, params);
		}
	}
}

} /* namespace rtabmap */
