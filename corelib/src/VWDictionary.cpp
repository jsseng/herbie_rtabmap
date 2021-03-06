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

#include "rtabmap/core/VWDictionary.h"
#include "rtabmap/core/VisualWord.h"

#include "rtabmap/core/Signature.h"
#include "rtabmap/core/DBDriver.h"
#include "rtabmap/core/Parameters.h"
#include "rtabmap/core/FlannIndex.h"

#include "rtabmap/utilite/UtiLite.h"

#include <opencv2/opencv_modules.hpp>
#include <chrono>

#if CV_MAJOR_VERSION < 3
#ifdef HAVE_OPENCV_GPU
#include <opencv2/gpu/gpu.hpp>
#endif
#else
#include <opencv2/core/cuda.hpp>
#ifdef HAVE_OPENCV_CUDAFEATURES2D
#include <opencv2/cudafeatures2d.hpp>
#endif
#endif

#include <fstream>
#include <string>

#define KDTREE_SIZE 4
#define KNN_CHECKS 32

namespace rtabmap
{

const int VWDictionary::ID_START = 1;
const int VWDictionary::ID_INVALID = 0;

VWDictionary::VWDictionary(const ParametersMap & parameters) :
	_totalActiveReferences(0),
	_incrementalDictionary(Parameters::defaultKpIncrementalDictionary()),
	_incrementalFlann(Parameters::defaultKpIncrementalFlann()),
	_rebalancingFactor(Parameters::defaultKpFlannRebalancingFactor()),
	_byteToFloat(Parameters::defaultKpByteToFloat()),
	_nndrRatio(Parameters::defaultKpNndrRatio()),
	_newDictionaryPath(Parameters::defaultKpDictionaryPath()),
	_newWordsComparedTogether(Parameters::defaultKpNewWordsComparedTogether()),
	_lastWordId(0),
	useDistanceL1_(false),
	_flannIndex(new FlannIndex()),
	_strategy(kNNBruteForce)
{
	this->setNNStrategy((NNStrategy)Parameters::defaultKpNNStrategy());
	this->parseParameters(parameters);
	this->build_cache = 0;
	this->read_cache = 0;
}

VWDictionary::~VWDictionary()
{
	this->clear();
	delete _flannIndex;
}

void VWDictionary::parseParameters(const ParametersMap & parameters)
{
	ParametersMap::const_iterator iter;
	Parameters::parse(parameters, Parameters::kKpNndrRatio(), _nndrRatio);
	Parameters::parse(parameters, Parameters::kKpNewWordsComparedTogether(), _newWordsComparedTogether);
	Parameters::parse(parameters, Parameters::kKpIncrementalFlann(), _incrementalFlann);
	Parameters::parse(parameters, Parameters::kKpFlannRebalancingFactor(), _rebalancingFactor);
	bool byteToFloat = _byteToFloat;
	Parameters::parse(parameters, Parameters::kKpByteToFloat(), _byteToFloat);

	UASSERT_MSG(_nndrRatio > 0.0f, uFormat("String=%s value=%f", uContains(parameters, Parameters::kKpNndrRatio())?parameters.at(Parameters::kKpNndrRatio()).c_str():"", _nndrRatio).c_str());

	bool incrementalDictionary = _incrementalDictionary;
	if((iter=parameters.find(Parameters::kKpDictionaryPath())) != parameters.end())
	{
		_newDictionaryPath = (*iter).second.c_str();
	}
	if((iter=parameters.find(Parameters::kKpIncrementalDictionary())) != parameters.end())
	{
		incrementalDictionary = uStr2Bool((*iter).second.c_str());
	}

	// Verifying hypotheses strategy
	bool treeUpdated = false;
	if((iter=parameters.find(Parameters::kKpNNStrategy())) != parameters.end())
	{
		NNStrategy nnStrategy = (NNStrategy)std::atoi((*iter).second.c_str());
		treeUpdated = this->setNNStrategy(nnStrategy);
	}
	if(!treeUpdated && byteToFloat!=_byteToFloat && _strategy == kNNFlannKdTree)
	{
		UINFO("KDTree: Binary to Float conversion approach has changed, re-initialize kd-tree.");
		_dataTree = cv::Mat();
		_notIndexedWords = uKeysSet(_visualWords);
		_removedIndexedWords.clear();
		this->update();
	}

	if(incrementalDictionary)
	{
		this->setIncrementalDictionary();
	}
	_incrementalDictionary = incrementalDictionary;
}

void VWDictionary::setIncrementalDictionary()
{
	if(!_incrementalDictionary)
	{
		_incrementalDictionary = true;
		if(_visualWords.size())
		{
			UWARN("Incremental dictionary set: already loaded visual words (%d) from the fixed dictionary will be included in the incremental one.", _visualWords.size());
		}
	}
	_dictionaryPath = "";
	_newDictionaryPath = "";
}

void VWDictionary::setFixedDictionary(const std::string & dictionaryPath)
{
	UDEBUG("");
	if(!dictionaryPath.empty())
	{
		if((!_incrementalDictionary && _dictionaryPath.compare(dictionaryPath) != 0) ||
		   _visualWords.size() == 0)
		{
			UINFO("incremental=%d, oldPath=%s newPath=%s, visual words=%d",
					_incrementalDictionary?1:0, _dictionaryPath.c_str(), dictionaryPath.c_str(), (int)_visualWords.size());

			if(UFile::getExtension(dictionaryPath).compare("db") == 0)
			{
				UWARN("Loading fixed vocabulary \"%s\", this may take a while...", dictionaryPath.c_str());
				DBDriver * driver = DBDriver::create();
				if(driver->openConnection(dictionaryPath, false))
				{
					driver->load(this, false);
					for(std::map<int, VisualWord*>::iterator iter=_visualWords.begin(); iter!=_visualWords.end(); ++iter)
					{
						iter->second->setSaved(true);
					}
					_incrementalDictionary = _visualWords.size()==0;
					driver->closeConnection(false);
				}
				else
				{
					UERROR("Could not load dictionary from database %s", dictionaryPath.c_str());
				}
				delete driver;
			}
			else
			{
				UWARN("Loading fixed vocabulary \"%s\", this may take a while...", dictionaryPath.c_str());
				std::ifstream file;
				file.open(dictionaryPath.c_str(), std::ifstream::in);
				if(file.good())
				{
					UDEBUG("Deleting old dictionary and loading the new one from \"%s\"", dictionaryPath.c_str());
					UTimer timer;

					// first line is the header
					std::string str;
					std::list<std::string> strList;
					std::getline(file, str);
					strList = uSplitNumChar(str);
					int dimension  = 0;
					for(std::list<std::string>::iterator iter = strList.begin(); iter != strList.end(); ++iter)
					{
						if(uIsDigit(iter->at(0)))
						{
							dimension = std::atoi(iter->c_str());
							break;
						}
					}
					UDEBUG("descriptor dimension = %d", dimension);

					if(dimension <= 0 || dimension > 1000)
					{
						UERROR("Invalid dictionary file, visual word dimension (%d) is not valid, \"%s\"", dimension, dictionaryPath.c_str());
					}
					else
					{
						// Process all words
						while(file.good())
						{
							std::getline(file, str);
							strList = uSplit(str);
							if((int)strList.size() == dimension+1)
							{
								//first one is the visual word id
								std::list<std::string>::iterator iter = strList.begin();
								int id = std::atoi(iter->c_str());
								cv::Mat descriptor(1, dimension, CV_32F);
								++iter;
								int i=0;

								//get descriptor
								for(;i<dimension && iter != strList.end(); ++i, ++iter)
								{
									descriptor.at<float>(i) = uStr2Float(*iter);
								}
								if(i != dimension)
								{
									UERROR("Loaded word has not the same size (%d) than descriptor size previously detected (%d).", i, dimension);
								}

								VisualWord * vw = new VisualWord(id, descriptor, 0);
								vw->setSaved(true);
								_visualWords.insert(_visualWords.end(), std::pair<int, VisualWord*>(id, vw));
								_notIndexedWords.insert(_notIndexedWords.end(), id);
								_unusedWords.insert(_unusedWords.end(), std::pair<int, VisualWord*>(id, vw));
							}
							else if(!str.empty())
							{
								UWARN("Cannot parse line \"%s\"", str.c_str());
							}
						}
						if(_visualWords.size())
						{
							UWARN("Loaded %d words!", (int)_visualWords.size());
						}
					}
				}
				else
				{
					UERROR("Cannot open dictionary file \"%s\"", dictionaryPath.c_str());
				}
				file.close();
			}

			if(_visualWords.size() == 0)
			{
				_incrementalDictionary = _visualWords.size()==0;
				UWARN("No words loaded, cannot set a fixed dictionary.", (int)_visualWords.size());
			}
			else
			{
				_dictionaryPath = dictionaryPath;
				_newDictionaryPath = dictionaryPath;
				_incrementalDictionary = false;
				this->update();
				UWARN("Loaded %d words!", (int)_visualWords.size());
				std::cout << "JS----VWDictionary.cpp - setFixedDictionary() - loaded words complete----";
			}
		}
		else if(!_incrementalDictionary)
		{
			UDEBUG("Dictionary \"%s\" already loaded...", dictionaryPath.c_str());
		}
		else
		{
			UERROR("Cannot change to a fixed dictionary if there are already words (%d) in the incremental one.", _visualWords.size());
		}
	}
	else if(_incrementalDictionary && _visualWords.size())
	{
		UWARN("Cannot change to fixed dictionary, %d words already loaded as incremental", (int)_visualWords.size());
	}
	else
	{
		_incrementalDictionary = false;
	}
	_dictionaryPath = dictionaryPath;
	_newDictionaryPath = dictionaryPath;
}

bool VWDictionary::setNNStrategy(NNStrategy strategy)
{
#if CV_MAJOR_VERSION < 3
#ifdef HAVE_OPENCV_GPU
	if(strategy == kNNBruteForceGPU && !cv::gpu::getCudaEnabledDeviceCount())
	{
		UERROR("Nearest neighobr strategy \"kNNBruteForceGPU\" chosen but no CUDA devices found! Doing \"kNNBruteForce\" instead.");
		strategy = kNNBruteForce;
	}
#else
	if(strategy == kNNBruteForceGPU)
	{
		UERROR("Nearest neighobr strategy \"kNNBruteForceGPU\" chosen but OpenCV is not built with GPU/cuda module! Doing \"kNNBruteForce\" instead.");
		strategy = kNNBruteForce;
	}
#endif
#else
#ifdef HAVE_OPENCV_CUDAFEATURES2D
	if(strategy == kNNBruteForceGPU && !cv::cuda::getCudaEnabledDeviceCount())
	{
		UERROR("Nearest neighobr strategy \"kNNBruteForceGPU\" chosen but no CUDA devices found! Doing \"kNNBruteForce\" instead.");
		strategy = kNNBruteForce;
	}
#else
	if(strategy == kNNBruteForceGPU)
	{
		UERROR("Nearest neighobr strategy \"kNNBruteForceGPU\" chosen but OpenCV cudafeatures2d module is not found! Doing \"kNNBruteForce\" instead.");
		strategy = kNNBruteForce;
	}
#endif
#endif

	if(strategy>=kNNUndef)
	{
		UERROR("Nearest neighobr strategy \"%d\" chosen but this strategy cannot be used with a dictionary! Doing \"kNNBruteForce\" instead.");
		strategy = kNNBruteForce;
	}

	bool update = _strategy != strategy;
	_strategy = strategy;
	if(update)
	{
		if(_notIndexedWords.size() != _visualWords.size() || !_dataTree.empty())
		{
			UINFO("Nearest neighbor strategy has changed, re-initialize search tree.");
		}
		_dataTree = cv::Mat();
		_notIndexedWords = uKeysSet(_visualWords);
		_removedIndexedWords.clear();
                std::cout << "---------Calling update 1--------------------" << std::endl;
		this->update();
		return true;
	}
	return false;
}

int VWDictionary::getLastIndexedWordId() const
{
	if(_mapIndexId.size())
	{
		return _mapIndexId.rbegin()->second;
	}
	else
	{
		return 0;
	}
}

unsigned int VWDictionary::getIndexedWordsCount() const
{
	return _flannIndex->indexedFeatures();
}

unsigned int VWDictionary::getIndexMemoryUsed() const
{
	return _flannIndex->memoryUsed();
}

unsigned long VWDictionary::getMemoryUsed() const
{
	long memoryUsage = sizeof(VWDictionary);
	memoryUsage += getIndexMemoryUsed();
	memoryUsage += _dataTree.total()*_dataTree.elemSize();
	if(!_visualWords.empty())
	{
		memoryUsage += _visualWords.size()*(sizeof(int) + _visualWords.rbegin()->second->getMemoryUsed() + sizeof(std::map<int, VisualWord *>::iterator)) + sizeof(std::map<int, VisualWord *>);
		if(_dataTree.empty() &&
			_visualWords.begin()->second->getDescriptor().type() == CV_8U &&
			_strategy == kNNFlannKdTree)
		{
			// Binary descriptors were converted to float, and not included in _dataTree
			memoryUsage += _visualWords.size() * _visualWords.begin()->second->getDescriptor().total() * sizeof(float) * (_byteToFloat?1:8);
		}
	}
	if(!_unusedWords.empty())
	{
		// they are the same words than in _visualWords, so just add the pointer size
		memoryUsage += _unusedWords.size()*(sizeof(int) + sizeof(VisualWord *)+sizeof(std::map<int, VisualWord *>::iterator)) + sizeof(std::map<int, VisualWord *>);
	}
	memoryUsage += _mapIndexId.size() * (sizeof(int)*2+sizeof(std::map<int ,int>::iterator)) + sizeof(std::map<int ,int>);
	memoryUsage += _mapIdIndex.size() * (sizeof(int)*2+sizeof(std::map<int ,int>::iterator)) + sizeof(std::map<int ,int>);
	memoryUsage += _notIndexedWords.size() * (sizeof(int)+sizeof(std::set<int>::iterator)) + sizeof(std::set<int>);
	memoryUsage += _removedIndexedWords.size() * (sizeof(int)+sizeof(std::set<int>::iterator)) + sizeof(std::set<int>);
	return memoryUsage;
}

cv::Mat VWDictionary::convertBinTo32F(const cv::Mat & descriptorsIn, bool byteToFloat)
{
	if(byteToFloat)
	{
		// Old approach
		cv::Mat descriptorsOut;
		descriptorsIn.convertTo(descriptorsOut, CV_32F);
		return descriptorsOut;
	}
	else
	{
		// New approach
		UASSERT(descriptorsIn.type() == CV_8UC1);
		cv::Mat descriptorsOut(descriptorsIn.rows, descriptorsIn.cols*8, CV_32FC1);
		for(int i=0; i<descriptorsIn.rows; ++i)
		{
			const unsigned char * ptrIn = descriptorsIn.ptr(i);
			float * ptrOut = descriptorsOut.ptr<float>(i);
			for(int j=0; j<descriptorsIn.cols; ++j)
			{
				int jo = j*8;
				ptrOut[jo] = (ptrIn[j] & 1) == 1?1.0f:0.0f;
				ptrOut[jo+1] = (ptrIn[j] & (1<<1)) != 0?1.0f:0.0f;
				ptrOut[jo+2] = (ptrIn[j] & (1<<2)) != 0?1.0f:0.0f;
				ptrOut[jo+3] = (ptrIn[j] & (1<<3)) != 0?1.0f:0.0f;
				ptrOut[jo+4] = (ptrIn[j] & (1<<4)) != 0?1.0f:0.0f;
				ptrOut[jo+5] = (ptrIn[j] & (1<<5)) != 0?1.0f:0.0f;
				ptrOut[jo+6] = (ptrIn[j] & (1<<6)) != 0?1.0f:0.0f;
				ptrOut[jo+7] = (ptrIn[j] & (1<<7)) != 0?1.0f:0.0f;
			}
		}
		return descriptorsOut;
	}
}

cv::Mat VWDictionary::convert32FToBin(const cv::Mat & descriptorsIn, bool byteToFloat)
{
	if(byteToFloat)
	{
		// Old approach
		cv::Mat descriptorsOut;
		descriptorsIn.convertTo(descriptorsOut, CV_8UC1);
		return descriptorsOut;
	}
	else
	{
		// New approach
		UASSERT(descriptorsIn.type() == CV_32FC1 && descriptorsIn.cols % 8 == 0);
		cv::Mat descriptorsOut(descriptorsIn.rows, descriptorsIn.cols/8, CV_8UC1);
		for(int i=0; i<descriptorsIn.rows; ++i)
		{
			const float * ptrIn = descriptorsIn.ptr<float>(i);
			unsigned char * ptrOut = descriptorsOut.ptr(i);
			for(int j=0; j<descriptorsOut.cols; ++j)
			{
				int jo = j*8;
				ptrOut[j] =
						(unsigned char)(ptrIn[jo] == 0?0:1) |
						(ptrIn[jo+1] == 0?0:(1<<1)) |
						(ptrIn[jo+2] == 0?0:(1<<2)) |
						(ptrIn[jo+3] == 0?0:(1<<3)) |
						(ptrIn[jo+4] == 0?0:(1<<4)) |
						(ptrIn[jo+5] == 0?0:(1<<5)) |
						(ptrIn[jo+6] == 0?0:(1<<6)) |
						(ptrIn[jo+7] == 0?0:(1<<7));
			}
		}
		return descriptorsOut;
	}
}

void VWDictionary::update()
{
	ULOGGER_DEBUG("incremental=%d", _incrementalDictionary?1:0);
	if(!_incrementalDictionary)
	{
		std::cout << "-----------------fixed dictionary--------------" << std::endl;
		// reload the fixed dictionary if it has been cleared or not yet initialized
		this->setFixedDictionary(_newDictionaryPath);

		if(!_incrementalDictionary && !_notIndexedWords.size())
		{
			// No need to update the search index if we
			// use a fixed dictionary and the index is
			// already built
			return;
		}
	}

	if(_notIndexedWords.size() || _visualWords.size() == 0 || _removedIndexedWords.size())
	{
		if(_incrementalFlann &&
		   _strategy < kNNBruteForce &&
		   _visualWords.size())
		{
			ULOGGER_DEBUG("Incremental FLANN: Removing %d words...", (int)_removedIndexedWords.size());
			for(std::set<int>::iterator iter=_removedIndexedWords.begin(); iter!=_removedIndexedWords.end(); ++iter)
			{
				UASSERT(uContains(_mapIdIndex, *iter));
				UASSERT(uContains(_mapIndexId, _mapIdIndex.at(*iter)));
				_flannIndex->removePoint(_mapIdIndex.at(*iter));
				_mapIndexId.erase(_mapIdIndex.at(*iter));
				_mapIdIndex.erase(*iter);
				std::cout << "---------removing word-------------" << std::endl;
			}
			ULOGGER_DEBUG("Incremental FLANN: Removing %d words... done!", (int)_removedIndexedWords.size());

			if(_notIndexedWords.size())
			{
				ULOGGER_DEBUG("Incremental FLANN: Inserting %d words...", (int)_notIndexedWords.size());
				for(std::set<int>::iterator iter=_notIndexedWords.begin(); iter!=_notIndexedWords.end(); ++iter)
				{
					VisualWord* w = uValue(_visualWords, *iter, (VisualWord*)0);
					UASSERT(w);

					cv::Mat descriptor;
					if(w->getDescriptor().type() == CV_8U)
					{
						useDistanceL1_ = true;
						if(_strategy == kNNFlannKdTree)
						{
							//in the flann index, data is stored as floats
							// unsigned char *d = (unsigned char*) (w->getDescriptor().data);
							// for (int i=0;i<256;i++) {
							// 	std::cout << w->getDescriptor().size();
							// 	d++;
							// }
							descriptor = convertBinTo32F(w->getDescriptor(), _byteToFloat);
						}
						else
						{
							descriptor = w->getDescriptor();
						}
					}
					else
					{
						descriptor = w->getDescriptor();
					}

					int index = 0;
					if(!_flannIndex->isBuilt())
					{
						UDEBUG("Building FLANN index...");
						switch(_strategy)
						{
						case kNNFlannNaive:
							_flannIndex->buildLinearIndex(descriptor, useDistanceL1_, _rebalancingFactor);
							break;
						case kNNFlannKdTree:
							UASSERT_MSG(descriptor.type() == CV_32F, "To use KdTree dictionary, float descriptors are required!");
							_flannIndex->buildKDTreeIndex(descriptor, KDTREE_SIZE, useDistanceL1_, _rebalancingFactor);
							break;
						case kNNFlannLSH:
							UASSERT_MSG(descriptor.type() == CV_8U, "To use LSH dictionary, binary descriptors are required!");
							_flannIndex->buildLSHIndex(descriptor, 12, 20, 2, _rebalancingFactor);
							break;
						default:
							UFATAL("Not supposed to be here!");
							break;
						}
						UDEBUG("Building FLANN index... done!");
					}
					else
					{
						UASSERT(descriptor.cols == _flannIndex->featuresDim());
						UASSERT(descriptor.type() == _flannIndex->featuresType());
						UASSERT(descriptor.rows == 1);
						index = _flannIndex->addPoints(descriptor).front();
					}
					std::pair<std::map<int, int>::iterator, bool> inserted;
					inserted = _mapIndexId.insert(std::pair<int, int>(index, w->id()));
					UASSERT(inserted.second);
					inserted = _mapIdIndex.insert(std::pair<int, int>(w->id(), index));
					UASSERT(inserted.second);
				}
				std::cout << "---(update() VWDictionary.cpp) - inserting words done!---" << std::endl;
				ULOGGER_DEBUG("Incremental FLANN: Inserting %d words... done!", (int)_notIndexedWords.size());
			}
		}
		else if(_strategy >= kNNBruteForce &&
				_notIndexedWords.size() &&
				_removedIndexedWords.size() == 0 &&
				_visualWords.size() &&
				_dataTree.rows)
		{
			//exit(0);
			std::cout << "---(update() VWDictionary.cpp) - else if 1---" << std::endl;
			//just add not indexed words
			int i = _dataTree.rows;
			_dataTree.reserve(_dataTree.rows + _notIndexedWords.size());
			for(std::set<int>::iterator iter=_notIndexedWords.begin(); iter!=_notIndexedWords.end(); ++iter)
			{
				VisualWord* w = uValue(_visualWords, *iter, (VisualWord*)0);
				UASSERT(w);
				UASSERT(w->getDescriptor().cols == _dataTree.cols);
				UASSERT(w->getDescriptor().type() == _dataTree.type());
				_dataTree.push_back(w->getDescriptor());
				_mapIndexId.insert(_mapIndexId.end(), std::pair<int, int>(i, w->id()));
				std::pair<std::map<int, int>::iterator, bool> inserted = _mapIdIndex.insert(std::pair<int, int>(w->id(), i));
				UASSERT(inserted.second);
				++i;
			}
		}
		else
		{
			//exit(0);
			std::cout << "---(update() VWDictionary.cpp) - else 2---" << std::endl;
			_mapIndexId.clear();
			_mapIdIndex.clear();
			_dataTree = cv::Mat();
			_flannIndex->release();

			if(_visualWords.size())
			{
				std::cout << "---(update() VWDictionary.cpp) - else 3---" << std::endl;
				UTimer timer;
				timer.start();

				int dim = _visualWords.begin()->second->getDescriptor().cols;
				int type;
				if(_visualWords.begin()->second->getDescriptor().type() == CV_8U)
				{
					useDistanceL1_ = true;
					if(_strategy == kNNFlannKdTree)
					{
						type = CV_32F;
						if(!_byteToFloat)
						{
							dim *= 8;
						}
					}
					else
					{
						type = _visualWords.begin()->second->getDescriptor().type();
					}
				}
				else
				{
					type = _visualWords.begin()->second->getDescriptor().type();
				}

				UASSERT(type == CV_32F || type == CV_8U);
				UASSERT(dim > 0);

				// Create the data matrix
				_dataTree = cv::Mat(_visualWords.size(), dim, type); // SURF descriptors are CV_32F
				std::map<int, VisualWord*>::const_iterator iter = _visualWords.begin();
				for(unsigned int i=0; i < _visualWords.size(); ++i, ++iter)
				{
					cv::Mat descriptor;
					if(iter->second->getDescriptor().type() == CV_8U)
					{
						if(_strategy == kNNFlannKdTree)
						{
							descriptor = convertBinTo32F(iter->second->getDescriptor(), _byteToFloat);
						}
						else
						{
							descriptor = iter->second->getDescriptor();
						}
					}
					else
					{
						descriptor = iter->second->getDescriptor();
					}

					UASSERT_MSG(descriptor.type() == type, uFormat("%d vs %d", descriptor.type(), type).c_str());
					UASSERT_MSG(descriptor.cols == dim, uFormat("%d vs %d", descriptor.cols, dim).c_str());

					descriptor.copyTo(_dataTree.row(i));
					_mapIndexId.insert(_mapIndexId.end(), std::pair<int, int>(i, iter->second->id()));
					_mapIdIndex.insert(_mapIdIndex.end(), std::pair<int, int>(iter->second->id(), i));
				}

				ULOGGER_DEBUG("_mapIndexId.size() = %d, words.size()=%d, _dim=%d",_mapIndexId.size(), _visualWords.size(), dim);
				ULOGGER_DEBUG("copying data = %f s", timer.ticks());

				switch(_strategy)
				{
				case kNNFlannNaive:
					_flannIndex->buildLinearIndex(_dataTree, useDistanceL1_, _incrementalDictionary&&_incrementalFlann?_rebalancingFactor:1);
					break;
				case kNNFlannKdTree:
					UASSERT_MSG(type == CV_32F, "To use KdTree dictionary, float descriptors are required!");
					_flannIndex->buildKDTreeIndex(_dataTree, KDTREE_SIZE, useDistanceL1_, _incrementalDictionary&&_incrementalFlann?_rebalancingFactor:1);
					break;
				case kNNFlannLSH:
					UASSERT_MSG(type == CV_8U, "To use LSH dictionary, binary descriptors are required!");
					_flannIndex->buildLSHIndex(_dataTree, 12, 20, 2, _incrementalDictionary&&_incrementalFlann?_rebalancingFactor:1);
					break;
				default:
					break;
				}

				ULOGGER_DEBUG("Time to create kd tree = %f s", timer.ticks());
			}
		}
		UDEBUG("Dictionary updated! (size=%d added=%d removed=%d)",
				_dataTree.rows, _notIndexedWords.size(), _removedIndexedWords.size());
	}
	else
	{
		UDEBUG("Dictionary has not changed, so no need to update it! (size=%d)", _dataTree.rows);
	}
	_notIndexedWords.clear();
	_removedIndexedWords.clear();
	UDEBUG("");
}

void VWDictionary::clear(bool printWarningsIfNotEmpty)
{
	ULOGGER_DEBUG("");
	if(printWarningsIfNotEmpty)
	{
		if(_visualWords.size() && _incrementalDictionary)
		{
			UWARN("Visual dictionary would be already empty here (%d words still in dictionary).", (int)_visualWords.size());
		}
		if(_notIndexedWords.size())
		{
			UWARN("Not indexed words should be empty here (%d words still not indexed)", (int)_notIndexedWords.size());
		}
	}
	for(std::map<int, VisualWord *>::iterator i=_visualWords.begin(); i!=_visualWords.end(); ++i)
	{
		delete (*i).second;
	}
	_visualWords.clear();
	_notIndexedWords.clear();
	_removedIndexedWords.clear();
	_totalActiveReferences = 0;
	_lastWordId = 0;
	_dataTree = cv::Mat();
	_mapIndexId.clear();
	_mapIdIndex.clear();
	_unusedWords.clear();
	_flannIndex->release();
	useDistanceL1_ = false;
}

int VWDictionary::getNextId()
{
	return ++_lastWordId;
}

void VWDictionary::addWordRef(int wordId, int signatureId)
{
	VisualWord * vw = 0;
	vw = uValue(_visualWords, wordId, vw);
	if(vw)
	{
		vw->addRef(signatureId);
		_totalActiveReferences += 1;

		_unusedWords.erase(vw->id());
	}
	else
	{
		UERROR("Not found word %d (dict size=%d)", wordId, (int)_visualWords.size());
	}
}

void VWDictionary::removeAllWordRef(int wordId, int signatureId)
{
	VisualWord * vw = 0;
	vw = uValue(_visualWords, wordId, vw);
	if(vw)
	{
		_totalActiveReferences -= vw->removeAllRef(signatureId);
		if(vw->getReferences().size() == 0)
		{
			_unusedWords.insert(std::pair<int, VisualWord*>(vw->id(), vw));
		}
	}
}

std::list<int> VWDictionary::addNewWords(
		const cv::Mat & descriptorsIn,
		int signatureId)
{
	UDEBUG("id=%d descriptors=%d", signatureId, descriptorsIn.rows);
	UTimer timer;
	std::list<int> wordIds;
	if(descriptorsIn.rows == 0 || descriptorsIn.cols == 0)
	{
		UERROR("Descriptors size is null!");
		return wordIds;
	}

	if(!_incrementalDictionary && _visualWords.empty())
	{
		UERROR("Dictionary mode is set to fixed but no words are in it!");
		return wordIds;
	}

	// verify we have the same features
	int dim = 0;
	int type = -1;
	if(_visualWords.size())
	{
		dim = _visualWords.begin()->second->getDescriptor().cols;
		type = _visualWords.begin()->second->getDescriptor().type();
		UASSERT(type == CV_32F || type == CV_8U);
	}
	if(dim && dim != descriptorsIn.cols)
	{
		UERROR("Descriptors (size=%d) are not the same size as already added words in dictionary(size=%d)", descriptorsIn.cols, dim);
		return wordIds;
	}
	if(type>=0 && type != descriptorsIn.type())
	{
		UERROR("Descriptors (type=%d) are not the same type as already added words in dictionary(type=%d)", descriptorsIn.type(), type);
		return wordIds;
	}

	// now compare with the actual index
	cv::Mat descriptors;
	if(descriptorsIn.type() == CV_8U)
	{
		useDistanceL1_ = true;
		if(_strategy == kNNFlannKdTree)
		{
			descriptors = convertBinTo32F(descriptorsIn, _byteToFloat);
		}
		else
		{
			descriptors = descriptorsIn;
		}
	}
	else
	{
		descriptors = descriptorsIn;
	}
	dim = 0;
	type = -1;
	if(_dataTree.rows || _flannIndex->isBuilt())
	{
		dim = _flannIndex->isBuilt()?_flannIndex->featuresDim():_dataTree.cols;
		type = _flannIndex->isBuilt()?_flannIndex->featuresType():_dataTree.type();
		UASSERT(type == CV_32F || type == CV_8U);
	}

	if(dim && dim != descriptors.cols)
	{
		UERROR("Descriptors (size=%d) are not the same size as already added words in dictionary(size=%d)", descriptors.cols, dim);
		return wordIds;
	}

	if(type>=0 && type != descriptors.type())
	{
		UERROR("Descriptors (type=%d) are not the same type as already added words in dictionary(type=%d)", descriptors.type(), type);
		return wordIds;
	}

	int dupWordsCountFromDict= 0;
	int dupWordsCountFromLast= 0;

	unsigned int k=2; // k nearest neighbors

	cv::Mat newWords;
	std::vector<int> newWordsId;

	cv::Mat results;
	cv::Mat dists;
	std::vector<std::vector<cv::DMatch> > matches;
	bool bruteForce = false;
	bool isL2NotSqr = false;

	UTimer timerLocal;
	timerLocal.start();

	if(_flannIndex->isBuilt() || (!_dataTree.empty() && _dataTree.rows >= (int)k))
	{
		//Find nearest neighbors
		UDEBUG("newPts.total()=%d ", descriptors.rows);

		if(_strategy == kNNFlannNaive || _strategy == kNNFlannKdTree || _strategy == kNNFlannLSH)
		{
			_flannIndex->knnSearch(descriptors, results, dists, k, KNN_CHECKS);
		}
		else if(_strategy == kNNBruteForce)
		{
			bruteForce = true;
			cv::BFMatcher matcher(descriptors.type()==CV_8U?cv::NORM_HAMMING:cv::NORM_L2SQR);
			matcher.knnMatch(descriptors, _dataTree, matches, k);
		}
		else if(_strategy == kNNBruteForceGPU)
		{
			bruteForce = true;
#if CV_MAJOR_VERSION < 3
#ifdef HAVE_OPENCV_GPU
			cv::gpu::GpuMat newDescriptorsGpu(descriptors);
			cv::gpu::GpuMat lastDescriptorsGpu(_dataTree);
			if(descriptors.type()==CV_8U)
			{
				cv::gpu::BruteForceMatcher_GPU<cv::Hamming> gpuMatcher;
				gpuMatcher.knnMatch(newDescriptorsGpu, lastDescriptorsGpu, matches, k);
			}
			else
			{
				cv::gpu::BruteForceMatcher_GPU<cv::L2<float> > gpuMatcher;
				gpuMatcher.knnMatch(newDescriptorsGpu, lastDescriptorsGpu, matches, k);
				isL2NotSqr = true;
			}
#else
			UERROR("Cannot use brute Force GPU because OpenCV is not built with gpu module.");
#endif
#else
#ifdef HAVE_OPENCV_CUDAFEATURES2D
			cv::cuda::GpuMat newDescriptorsGpu(descriptors);
			cv::cuda::GpuMat lastDescriptorsGpu(_dataTree);
			cv::Ptr<cv::cuda::DescriptorMatcher> gpuMatcher;
			if(descriptors.type()==CV_8U)
			{
				gpuMatcher = cv::cuda::DescriptorMatcher::createBFMatcher(cv::NORM_HAMMING);
				gpuMatcher->knnMatch(newDescriptorsGpu, lastDescriptorsGpu, matches, k);
			}
			else
			{
				gpuMatcher = cv::cuda::DescriptorMatcher::createBFMatcher(cv::NORM_L2);
				gpuMatcher->knnMatch(newDescriptorsGpu, lastDescriptorsGpu, matches, k);
				isL2NotSqr = true;
			}
#else
			UERROR("Cannot use brute Force GPU because OpenCV is not built with cuda module.");
#endif
#endif
		}
		else
		{
			UFATAL("");
		}

		// In case of binary descriptors
		if(dists.type() == CV_32S)
		{
			cv::Mat temp;
			dists.convertTo(temp, CV_32F);
			dists = temp;
		}

		UDEBUG("Time to find nn = %f s", timerLocal.ticks());
	}

	// Process results
	for(int i = 0; i < descriptors.rows; ++i)
	{
		std::multimap<float, int> fullResults; // Contains results from the kd-tree search and the naive search in new words
		if(!bruteForce && dists.cols)
		{
			for(int j=0; j<dists.cols; ++j)
			{
				float d = dists.at<float>(i,j);
				int index;
				if (sizeof(size_t) == 8)
				{
					index = *((size_t*)&results.at<double>(i, j));
				}
				else
				{
					index = *((size_t*)&results.at<int>(i, j));
				}
				int id = uValue(_mapIndexId, index);
				// int id = _mapIndexId[index];
				if(d >= 0.0f && id != 0)
				{
					fullResults.insert(std::pair<float, int>(d, id));
				}
				else
				{
					break;
				}
			}
		}
		else if(bruteForce && matches.size())
		{
			for(unsigned int j=0; j<matches.at(i).size(); ++j)
			{
				float d = matches.at(i).at(j).distance;
				int id = uValue(_mapIndexId, matches.at(i).at(j).trainIdx);
				// int id = _mapIndexId[matches.at(i).at(j).trainIdx];
				if(d >= 0.0f && id != 0)
				{
					if(isL2NotSqr)
					{
						// Make it compatible with L2SQR format of flann
						d*=d;
					}
					fullResults.insert(std::pair<float, int>(d, id));
				}
				else
				{
					break;
				}
			}
		}

		// Check if this descriptor matches with a word from the last signature (a word not already added to the tree)
		if(_newWordsComparedTogether && newWords.rows)
		{
			std::vector<std::vector<cv::DMatch> > matchesNewWords;
			cv::BFMatcher matcher(descriptors.type()==CV_8U?cv::NORM_HAMMING:useDistanceL1_?cv::NORM_L1:cv::NORM_L2SQR);
			UASSERT(descriptors.cols == newWords.cols && descriptors.type() == newWords.type());
			matcher.knnMatch(descriptors.row(i), newWords, matchesNewWords, newWords.rows>1?2:1);
			UASSERT(matchesNewWords.size() == 1);
			for(unsigned int j=0; j<matchesNewWords.at(0).size(); ++j)
			{
				float d = matchesNewWords.at(0).at(j).distance;
				int id = newWordsId[matchesNewWords.at(0).at(j).trainIdx];
				if(d >= 0.0f && id != 0)
				{
					fullResults.insert(std::pair<float, int>(d, id));
				}
				else
				{
					break;
				}
			}
		}

		if(_incrementalDictionary)
		{
			bool badDist = false;
			if(fullResults.size() == 0)
			{
				badDist = true;
			}
			if(!badDist)
			{
				if(fullResults.size() >= 2)
				{
					// Apply NNDR
					if(fullResults.begin()->first > _nndrRatio * (++fullResults.begin())->first)
					{
						badDist = true; // Rejected
					}
				}
				else
				{
					badDist = true; // Rejected
				}
			}

			if(badDist)
			{
				// use original descriptor
				VisualWord * vw = new VisualWord(getNextId(), descriptorsIn.row(i), signatureId);
				_visualWords.insert(_visualWords.end(), std::pair<int, VisualWord *>(vw->id(), vw));
				_notIndexedWords.insert(_notIndexedWords.end(), vw->id());
				newWords.push_back(descriptors.row(i));
				newWordsId.push_back(vw->id());
				wordIds.push_back(vw->id());
				UASSERT(vw->id()>0);
			}
			else
			{
				if(_notIndexedWords.find(fullResults.begin()->second) != _notIndexedWords.end())
				{
					++dupWordsCountFromLast;
				}
				else
				{
					++dupWordsCountFromDict;
				}

				this->addWordRef(fullResults.begin()->second, signatureId);
				wordIds.push_back(fullResults.begin()->second);
			}
		}
		else if(fullResults.size())
		{
			// If the dictionary is not incremental, just take the nearest word
			++dupWordsCountFromDict;
			this->addWordRef(fullResults.begin()->second, signatureId);
			wordIds.push_back(fullResults.begin()->second);
			UASSERT(fullResults.begin()->second>0);
		}
	}
	ULOGGER_DEBUG("naive search and add ref/words time = %f s", timerLocal.ticks());

	ULOGGER_DEBUG("%d new words added...", _notIndexedWords.size());
	ULOGGER_DEBUG("%d duplicated words added (from current image = %d)...",
			dupWordsCountFromDict+dupWordsCountFromLast, dupWordsCountFromLast);
	UDEBUG("total time %fs", timer.ticks());

	_totalActiveReferences += _notIndexedWords.size();
	return wordIds;
}

std::vector<int> VWDictionary::findNN(const std::list<VisualWord *> & vws) const
{
	UTimer timer;
	timer.start();

	if(_visualWords.size() && vws.size())
	{
		int type = (*vws.begin())->getDescriptor().type();
		int dim = (*vws.begin())->getDescriptor().cols;

		if(dim != _visualWords.begin()->second->getDescriptor().cols)
		{
			UERROR("Descriptors (size=%d) are not the same size as already added words in dictionary(size=%d)", (*vws.begin())->getDescriptor().cols, dim);
			return std::vector<int>(vws.size(), 0);
		}

		if(type != _visualWords.begin()->second->getDescriptor().type())
		{
			UERROR("Descriptors (type=%d) are not the same type as already added words in dictionary(type=%d)", (*vws.begin())->getDescriptor().type(), type);
			return std::vector<int>(vws.size(), 0);
		}

		// fill the request matrix
		int index = 0;
		VisualWord * vw;
		cv::Mat query(vws.size(), dim, type);
		for(std::list<VisualWord *>::const_iterator iter=vws.begin(); iter!=vws.end(); ++iter, ++index)
		{
			vw = *iter;
			UASSERT(vw);

			UASSERT(vw->getDescriptor().cols == dim);
			UASSERT(vw->getDescriptor().type() == type);

			vw->getDescriptor().copyTo(query.row(index));
		}
		ULOGGER_DEBUG("Preparation time = %fs", timer.ticks());

		return findNN(query);
	}
	return std::vector<int>(vws.size(), 0);
}

std::vector<int> VWDictionary::findNN(const cv::Mat & queryIn) const
{
	UTimer timer;
	timer.start();
	std::vector<int> resultIds(queryIn.rows, 0);
	unsigned int k=2; // k nearest neighbor

	if(_visualWords.size() && queryIn.rows)
	{
		// verify we have the same features
		int dim = _visualWords.begin()->second->getDescriptor().cols;
		int type = _visualWords.begin()->second->getDescriptor().type();
		UASSERT(type == CV_32F || type == CV_8U);

		if(dim != queryIn.cols)
		{
			UERROR("Descriptors (size=%d) are not the same size as already added words in dictionary(size=%d)", queryIn.cols, dim);
			return resultIds;
		}
		if(type != queryIn.type())
		{
			UERROR("Descriptors (type=%d) are not the same type as already added words in dictionary(type=%d)", queryIn.type(), type);
			return resultIds;
		}

		// now compare with the actual index
		cv::Mat query;
		if(queryIn.type() == CV_8U)
		{
			if(_strategy == kNNFlannKdTree)
			{
				query = convertBinTo32F(queryIn, _byteToFloat);
			}
			else
			{
				query = queryIn;
			}
		}
		else
		{
			query = queryIn;
		}
		dim = 0;
		type = -1;
		if(_dataTree.rows || _flannIndex->isBuilt())
		{
			dim = _flannIndex->isBuilt()?_flannIndex->featuresDim():_dataTree.cols;
			type = _flannIndex->isBuilt()?_flannIndex->featuresType():_dataTree.type();
			UASSERT(type == CV_32F || type == CV_8U);
		}

		if(dim && dim != query.cols)
		{
			UERROR("Descriptors (size=%d) are not the same size as already added words in dictionary(size=%d)", query.cols, dim);
			return resultIds;
		}

		if(type>=0 && type != query.type())
		{
			UERROR("Descriptors (type=%d) are not the same type as already added words in dictionary(type=%d)", query.type(), type);
			return resultIds;
		}

		std::vector<std::vector<cv::DMatch> > matches;
		bool bruteForce = false;
		bool isL2NotSqr = false;
		cv::Mat results;
		cv::Mat dists;

		if(_flannIndex->isBuilt() || (!_dataTree.empty() && _dataTree.rows >= (int)k))
		{
			//Find nearest neighbors
			UDEBUG("query.rows=%d ", query.rows);

			if(_strategy == kNNFlannNaive || _strategy == kNNFlannKdTree || _strategy == kNNFlannLSH)
			{
				_flannIndex->knnSearch(query, results, dists, k, KNN_CHECKS);
			}
			else if(_strategy == kNNBruteForce)
			{
				bruteForce = true;
				cv::BFMatcher matcher(query.type()==CV_8U?cv::NORM_HAMMING:cv::NORM_L2SQR);
				matcher.knnMatch(query, _dataTree, matches, k);
			}
			else if(_strategy == kNNBruteForceGPU)
			{
				bruteForce = true;
#if CV_MAJOR_VERSION < 3
#ifdef HAVE_OPENCV_GPU
				cv::gpu::GpuMat newDescriptorsGpu(query);
				cv::gpu::GpuMat lastDescriptorsGpu(_dataTree);
				if(query.type()==CV_8U)
				{
					cv::gpu::BruteForceMatcher_GPU<cv::Hamming> gpuMatcher;
					gpuMatcher.knnMatch(newDescriptorsGpu, lastDescriptorsGpu, matches, k);
				}
				else
				{
					cv::gpu::BruteForceMatcher_GPU<cv::L2<float> > gpuMatcher;
					gpuMatcher.knnMatch(newDescriptorsGpu, lastDescriptorsGpu, matches, k);
					isL2NotSqr = true;
				}
#else
			UERROR("Cannot use brute Force GPU because OpenCV is not built with gpu module.");
#endif
#else
#ifdef HAVE_OPENCV_CUDAFEATURES2D
				cv::cuda::GpuMat newDescriptorsGpu(query);
				cv::cuda::GpuMat lastDescriptorsGpu(_dataTree);
				cv::cuda::GpuMat matchesGpu;
				cv::Ptr<cv::cuda::DescriptorMatcher> gpuMatcher;
				if(query.type()==CV_8U)
				{
					gpuMatcher = cv::cuda::DescriptorMatcher::createBFMatcher(cv::NORM_HAMMING);
					gpuMatcher->knnMatchAsync(newDescriptorsGpu, lastDescriptorsGpu, matchesGpu, k);
				}
				else
				{
					gpuMatcher = cv::cuda::DescriptorMatcher::createBFMatcher(cv::NORM_L2);
					gpuMatcher->knnMatchAsync(newDescriptorsGpu, lastDescriptorsGpu, matchesGpu, k);
					isL2NotSqr = true;
				}
				gpuMatcher->knnMatchConvert(matchesGpu, matches);
#else
			UERROR("Cannot use brute Force GPU because OpenCV is not built with cuda module.");
#endif
#endif
			}
			else
			{
				UFATAL("");
			}

			// In case of binary descriptors
			if(dists.type() == CV_32S)
			{
				cv::Mat temp;
				dists.convertTo(temp, CV_32F);
				dists = temp;
			}
		}
		ULOGGER_DEBUG("Search dictionary time = %fs", timer.ticks());

		std::map<int, int> mapIndexIdNotIndexed;
		std::vector<std::vector<cv::DMatch> > matchesNotIndexed;
		if(!_notIndexedWords.empty())
		{
			cv::Mat dataNotIndexed = cv::Mat::zeros(_notIndexedWords.size(), query.cols, query.type());
			unsigned int index = 0;
			VisualWord * vw;
			for(std::set<int>::iterator iter = _notIndexedWords.begin(); iter != _notIndexedWords.end(); ++iter, ++index)
			{
				vw = _visualWords.at(*iter);

				cv::Mat descriptor;
				if(vw->getDescriptor().type() == CV_8U)
				{
					if(_strategy == kNNFlannKdTree)
					{
						descriptor = convertBinTo32F(vw->getDescriptor(), _byteToFloat);
					}
					else
					{
						descriptor = vw->getDescriptor();
					}
				}
				else
				{
					descriptor = vw->getDescriptor();
				}
				UASSERT(vw != 0 && descriptor.cols == query.cols && descriptor.type() == query.type());
				descriptor.copyTo(dataNotIndexed.row(index));
				mapIndexIdNotIndexed.insert(mapIndexIdNotIndexed.end(), std::pair<int,int>(index, vw->id()));
			}
			// Find nearest neighbor
			ULOGGER_DEBUG("Searching in words not indexed...");
			cv::BFMatcher matcher(query.type()==CV_8U?cv::NORM_HAMMING:useDistanceL1_?cv::NORM_L1:cv::NORM_L2SQR);
			matcher.knnMatch(query, dataNotIndexed, matchesNotIndexed, dataNotIndexed.rows>1?2:1);
		}
		ULOGGER_DEBUG("Search not yet indexed words time = %fs", timer.ticks());

		for(int i=0; i<query.rows; ++i)
		{
			std::multimap<float, int> fullResults; // Contains results from the kd-tree search [and the naive search in new words]
			if(!bruteForce && dists.cols)
			{
				for(int j=0; j<dists.cols; ++j)
				{
					float d = dists.at<float>(i,j);
					int index;

					if (sizeof(size_t) == 8)
					{
						index = *((size_t*)&results.at<double>(i, j));
					}
					else
					{
						index = *((size_t*)&results.at<int>(i, j));
					}
					int id = uValue(_mapIndexId, index);
					// int id = _mapIndexId.at(index);
					if(d >= 0.0f && id != 0)
					{
						fullResults.insert(std::pair<float, int>(d, id));
					}
				}
			}
			else if(bruteForce && matches.size())
			{
				for(unsigned int j=0; j<matches.at(i).size(); ++j)
				{
					float d = matches.at(i).at(j).distance;
					int id = uValue(_mapIndexId, matches.at(i).at(j).trainIdx);
					// int id = _mapIndexId.at(matches.at(i).at(j).trainIdx);
					if(d >= 0.0f && id != 0)
					{
						if(isL2NotSqr)
						{
							// make it compatible with L2SQR from FLANN
							d*=d;
						}
						fullResults.insert(std::pair<float, int>(d, id));
					}
				}
			}

			// not indexed..
			if(matchesNotIndexed.size())
			{
				for(unsigned int j=0; j<matchesNotIndexed.at(i).size(); ++j)
				{
					float d = matchesNotIndexed.at(i).at(j).distance;
					int id = uValue(mapIndexIdNotIndexed, matchesNotIndexed.at(i).at(j).trainIdx);
					if(d >= 0.0f && id != 0)
					{
						fullResults.insert(std::pair<float, int>(d, id));
					}
					else
					{
						break;
					}
				}
			}

			if(_incrementalDictionary)
			{
				bool badDist = false;
				if(fullResults.size() == 0)
				{
					badDist = true;
				}
				if(!badDist)
				{
					if(fullResults.size() >= 2)
					{
						// Apply NNDR
						if(fullResults.begin()->first > _nndrRatio * (++fullResults.begin())->first)
						{
							badDist = true; // Rejected
						}
					}
					else
					{
						badDist = true; // Rejected
					}
				}

				if(!badDist)
				{
					resultIds[i] = fullResults.begin()->second; // Accepted
				}
			}
			else if(fullResults.size())
			{
				//Just take the nearest if the dictionary is not incremental
				resultIds[i] = fullResults.begin()->second; // Accepted
			}
		}
		ULOGGER_DEBUG("badDist check time = %fs", timer.ticks());
	}
	return resultIds;
}

void VWDictionary::addWord(VisualWord * vw)
{
	if(vw)
	{
		_visualWords.insert(std::pair<int, VisualWord *>(vw->id(), vw));
		_notIndexedWords.insert(vw->id());
		if(vw->getReferences().size())
		{
			_totalActiveReferences += uSum(uValues(vw->getReferences()));
		}
		else
		{
			_unusedWords.insert(std::pair<int, VisualWord *>(vw->id(), vw));
		}
		if(_lastWordId < vw->id())
		{
			_lastWordId = vw->id();
		}
	}
}

const VisualWord * VWDictionary::getWord(int id) const
{
	return uValue(_visualWords, id, (VisualWord *)0);
}

VisualWord * VWDictionary::getUnusedWord(int id) const
{
	return uValue(_unusedWords, id, (VisualWord *)0);
}

std::vector<VisualWord*> VWDictionary::getUnusedWords() const
{
	return uValues(_unusedWords);
}

std::vector<int> VWDictionary::getUnusedWordIds() const
{
	return uKeys(_unusedWords);
}

void VWDictionary::removeWords(const std::vector<VisualWord*> & words)
{
	//UDEBUG("Removing %d words from dictionary (current size=%d)", (int)words.size(), (int)_visualWords.size());
	for(unsigned int i=0; i<words.size(); ++i)
	{
		_visualWords.erase(words[i]->id());
		_unusedWords.erase(words[i]->id());
		if(_notIndexedWords.erase(words[i]->id()) == 0)
		{
			_removedIndexedWords.insert(words[i]->id());
		}
	}
}

void VWDictionary::deleteUnusedWords()
{
	std::vector<VisualWord*> unusedWords = uValues(_unusedWords);
	removeWords(unusedWords);
	for(unsigned int i=0; i<unusedWords.size(); ++i)
	{
		delete unusedWords[i];
	}
}

void VWDictionary::exportDictionary(const char * fileNameReferences, const char * fileNameDescriptors) const
{
	UDEBUG("");
	if(_visualWords.empty())
	{
		UWARN("Dictionary is empty, cannot export it!");
		return;
	}
	if(_visualWords.begin()->second->getDescriptor().type() != CV_32FC1)
	{
		UERROR("Exporting binary descriptors is not implemented!");
		return;
	}
    FILE* foutRef = 0;
    FILE* foutDesc = 0;
#ifdef _MSC_VER
    fopen_s(&foutRef, fileNameReferences, "w");
    fopen_s(&foutDesc, fileNameDescriptors, "w");
#else
    foutRef = fopen(fileNameReferences, "w");
    foutDesc = fopen(fileNameDescriptors, "w");
#endif

    if(foutRef)
    {
    	fprintf(foutRef, "WordID SignaturesID...\n");
    }
	if(foutDesc)
	{
		if(_visualWords.begin() == _visualWords.end())
		{
			fprintf(foutDesc, "WordID Descriptors...\n");
		}
		else
		{
			UDEBUG("");
			fprintf(foutDesc, "WordID Descriptors...%d\n", (*_visualWords.begin()).second->getDescriptor().cols);
		}
	}

	UDEBUG("Export %d words...", _visualWords.size());
    for(std::map<int, VisualWord *>::const_iterator iter=_visualWords.begin(); iter!=_visualWords.end(); ++iter)
    {
    	// References
    	if(foutRef)
    	{
			fprintf(foutRef, "%d ", (*iter).first);
			const std::map<int, int> ref = (*iter).second->getReferences();
			for(std::map<int, int>::const_iterator jter=ref.begin(); jter!=ref.end(); ++jter)
			{
				for(int i=0; i<(*jter).second; ++i)
				{
					fprintf(foutRef, "%d ", (*jter).first);
				}
			}
			fprintf(foutRef, "\n");
    	}

    	//Descriptors
    	if(foutDesc)
    	{
			fprintf(foutDesc, "%d ", (*iter).first);
			const float * desc = (const float *)(*iter).second->getDescriptor().data;
			int dim = (*iter).second->getDescriptor().cols;

			for(int i=0; i<dim; i++)
			{
				fprintf(foutDesc, "%f ", desc[i]);
			}
			fprintf(foutDesc, "\n");
    	}
    }

	if(foutRef)
		fclose(foutRef);
	if(foutDesc)
		fclose(foutDesc);
}

void VWDictionary::debug()
{
	std::cout << _flannIndex << std::endl;
	std::cout << "num flann features dimensions: " << _flannIndex->featuresDim() << std::endl;
	std::cout << "flann bytes used: " << _flannIndex->memoryUsed() << std::endl;
	std::cout << "num indexed features: " << _flannIndex->indexedFeatures() << std::endl;
	std::cout << "VWDictionary total active references: " << getTotalActiveReferences() << std::endl;
	// std::cout << "_vwd: " << test1 << std::endl;

	_flannIndex->debug();
	std::cout << "_mapIndexId: " << _mapIndexId.size() << std::endl;
	std::cout << "_mapIndexId: " << _mapIdIndex.size() << std::endl;
	std::cout << "_unusedWords: " << _unusedWords.size() << std::endl;
	std::cout << "_notIndexedWords: " << _notIndexedWords.size() << std::endl;
	std::cout << "_removedIndexedWords: " << _removedIndexedWords.size() << std::endl;

	// std::cout << "datatree rows: " << _dataTree.rows << std::endl;
	// std::cout << "datatree cols: " << _dataTree.cols << std::endl;
}

void VWDictionary::check_vwdictionary() {}

void VWDictionary::save_vwdictionary()
{
	std::ofstream* outfile;
	std::ifstream infile;
	outfile = new std::ofstream();
	outfile->open("/home/jseng/vwdictionary.dat", std::ios::out | std::ios::binary | std::ios::trunc );

	int visualword_num = _visualWords.size();
	outfile->write(reinterpret_cast<char *>(&visualword_num), sizeof(int));
	for(std::map<int, VisualWord*>::iterator iter=_visualWords.begin(); iter!=_visualWords.end(); ++iter)
	{
		visualword_num = iter->first;
		outfile->write(reinterpret_cast<char*> (&visualword_num),sizeof(int));  //write out the id
		(iter->second)->save_visualword(outfile);
	}
	
	std::cout << "saved _incrementalDictionary: " << _incrementalDictionary << std::endl;
	outfile->write(reinterpret_cast<const char *> (&_incrementalDictionary),sizeof(bool)); //_incrementalDictionary

	std::cout << "saved _incrementalFlann: " << _incrementalFlann << std::endl;
	outfile->write(reinterpret_cast<const char *> (&_incrementalFlann),sizeof(bool)); //_incrementalFlann

	std::cout << "saved rebalancing factor: " << _rebalancingFactor << std::endl;
	outfile->write(reinterpret_cast<const char *> (&_rebalancingFactor),sizeof(float)); //_rebalancingFactor

	std::cout << "saved _byteToFloat: " << _byteToFloat << std::endl;
	outfile->write(reinterpret_cast<const char *> (&_byteToFloat),sizeof(bool)); //_byteToFloat

	std::cout << "saved _nndrRatio: " << _nndrRatio << std::endl;
	outfile->write(reinterpret_cast<const char *> (&_nndrRatio),sizeof(float)); //_nndrRatio

	int dictionarypath_length = _dictionaryPath.length();
	outfile->write(reinterpret_cast<const char *> (&dictionarypath_length),sizeof(int)); //length of the dictionary path
	outfile->write(reinterpret_cast<const char *> (_dictionaryPath.data()),dictionarypath_length); //_dictionaryPath

	int newdictionarypath_length = _newDictionaryPath.length();
	outfile->write(reinterpret_cast<const char *> (&newdictionarypath_length),sizeof(int)); //length of the new dictionary path
	outfile->write(reinterpret_cast<const char *> (_newDictionaryPath.data()),newdictionarypath_length); //_newDictionaryPath

	outfile->write(reinterpret_cast<const char *> (&_newWordsComparedTogether),sizeof(bool)); //_newWordsCompareTogether
	std::cout << "saved _newWordsComparedTogether: " << _newWordsComparedTogether << std::endl;

	outfile->write(reinterpret_cast<const char *> (&_lastWordId),sizeof(int)); //_lastWordId
	std::cout << "saved _lastWordId: " << _lastWordId << std::endl;

	outfile->write(reinterpret_cast<const char *> (&useDistanceL1_),sizeof(bool)); //useDistanceL1_
	std::cout << "saved useDistanceL1_: " << useDistanceL1_ << std::endl;

	int strat = _strategy;
	outfile->write(reinterpret_cast<const char *> (&strat),sizeof(int)); //useDistanceL1_
	std::cout << "saved _strategy: " << _strategy << std::endl;

	//std::map<int, int> _mapIndexId;
	int mapIndexId_temp = _mapIndexId.size();
	outfile->write(reinterpret_cast<char *>(&mapIndexId_temp), sizeof(int));
	for(std::map<int, int>::iterator iter=_mapIndexId.begin(); iter!=_mapIndexId.end(); ++iter)
	{
		mapIndexId_temp = iter->first;
		outfile->write(reinterpret_cast<char*> (&mapIndexId_temp),sizeof(int));
		mapIndexId_temp = iter->second;
		outfile->write(reinterpret_cast<char*> (&mapIndexId_temp),sizeof(int));
	}

	// std::map<int, int> _mapIdIndex;
	int mapIdIndex_temp = _mapIdIndex.size();
	outfile->write(reinterpret_cast<char *>(&mapIdIndex_temp), sizeof(int));
	for(std::map<int, int>::iterator iter=_mapIdIndex.begin(); iter!=_mapIdIndex.end(); ++iter)
	{
		mapIdIndex_temp = iter->first;
		outfile->write(reinterpret_cast<char*> (&mapIdIndex_temp),sizeof(int));
		mapIdIndex_temp = iter->second;
		outfile->write(reinterpret_cast<char*> (&mapIdIndex_temp),sizeof(int));
	}

	// std::map<int, VisualWord *> _unusedWords; //<id,VisualWord*>, note that these words stay in _visualWords
	int unused_visualword_num = _unusedWords.size();
	std::cout << "saving unused visual word size: " << unused_visualword_num << std::endl;
	outfile->write(reinterpret_cast<char *>(&unused_visualword_num), sizeof(int));
	for(std::map<int, VisualWord*>::iterator iter=_unusedWords.begin(); iter!=_unusedWords.end(); ++iter)
	{
		unused_visualword_num = iter->first;
		outfile->write(reinterpret_cast<char*> (&unused_visualword_num),sizeof(int));  //write out the id
		(iter->second)->save_visualword(outfile);
	}

	// std::set<int> _notIndexedWords;			  // Words that are not indexed in the dictionary
	int notIndexedWords_size = _notIndexedWords.size();
	std::cout << "saving not indexed visual word size: " << notIndexedWords_size << std::endl;
	int notindexed_temp;
	outfile->write(reinterpret_cast<char *>(&notIndexedWords_size), sizeof(int));
	for(std::set<int>::iterator iter=_notIndexedWords.begin(); iter!=_notIndexedWords.end(); ++iter)
	{
		notindexed_temp = (*iter);
		outfile->write(reinterpret_cast<char*> (&notindexed_temp),sizeof(int));
	}

	// std::set<int> _removedIndexedWords;		  // Words not anymore in the dictionary but still indexed in the dictionary
	int removedIndexedWords_size = _removedIndexedWords.size();
	std::cout << "saving removedIndexedWords size: " << removedIndexedWords_size << std::endl;
	int removed_temp;
	outfile->write(reinterpret_cast<char *>(&removedIndexedWords_size), sizeof(int));
	for(std::set<int>::iterator iter=_removedIndexedWords.begin(); iter!=_removedIndexedWords.end(); ++iter)
	{
		removed_temp = (*iter);
		outfile->write(reinterpret_cast<char*> (&removed_temp),sizeof(int));
	}

	outfile->close();

	_flannIndex->save_index();

	//check data
	// infile.open("vwdictionary.dat", std::ios::in | std::ios::binary);
	// infile.read(reinterpret_cast<char*> (&visualword_num),sizeof(int));
	// std::cout << "test data: " << visualword_num << std::endl;
	// infile.close();
}

void VWDictionary::load_vwdictionary() {
	auto t1 = std::chrono::high_resolution_clock::now();
	std::ifstream* infile;
	infile = new std::ifstream();
	infile->open("/home/jseng/vwdictionary.dat", std::ios::in | std::ios::binary);

	_visualWords.clear();
	int visualword_size;
	VisualWord* v;
	int id;
	infile->read(reinterpret_cast<char *>(&visualword_size), sizeof(int));  //read in the number of visual words
	// std::cout << "restored visualword_size: " << visualword_size << std::endl;
	for(int i=0; i<visualword_size; i++) 
	{
		cv::Mat d;	
		d = cv::Mat(1, 32, CV_8U);	
		infile->read(reinterpret_cast<char *>(&id), sizeof(int));  //read in the id
		v = new VisualWord(id,d);
		v->load_visualword(infile);
		_visualWords.insert(std::pair<int,VisualWord*>(id,v));
	}

	auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
	std::cout << "vwdictionary load time 1 in milliseconds: " << fp_ms.count() << std::endl;
	
	infile->read(reinterpret_cast<char *> (&_incrementalDictionary),sizeof(bool)); //_incrementalDictionary
	// std::cout << "restored _incrementalDictionary factor: " << _incrementalDictionary << std::endl;

	infile->read(reinterpret_cast<char *> (&_incrementalFlann),sizeof(bool)); //_incrementalFlann
	// std::cout << "restored _incrementalFlann: " << _incrementalFlann << std::endl;

	infile->read(reinterpret_cast<char *> (&_rebalancingFactor),sizeof(float)); //_rebalancingFactor
	// std::cout << "restored rebalancing factor: " << _rebalancingFactor << std::endl;

	infile->read(reinterpret_cast<char *> (&_byteToFloat),sizeof(bool)); //_byteToFloat
	// std::cout << "restored _byteToFloat: " << _byteToFloat << std::endl;

	infile->read(reinterpret_cast<char *> (&_nndrRatio),sizeof(float)); //_nndrRatio
	// std::cout << "restored _nndrRatio : " << _nndrRatio << std::endl;

	char temp_string[100];
	int dictionarypath_length;
	infile->read(reinterpret_cast<char *> (&dictionarypath_length),sizeof(int)); //length of the dictionary path
	infile->read(temp_string,dictionarypath_length); //_dictionaryPath
	temp_string[dictionarypath_length] = 0;
	_dictionaryPath = temp_string;

	int newdictionarypath_length;
	infile->read(reinterpret_cast<char *> (&newdictionarypath_length),sizeof(int)); //length of the new dictionary path
	infile->read(temp_string,newdictionarypath_length); //_newDictionaryPath
	temp_string[newdictionarypath_length] = 0;
	_newDictionaryPath = temp_string;

	infile->read(reinterpret_cast<char *> (&_newWordsComparedTogether),sizeof(bool)); //_newWordsCompareTogether
	// std::cout << "restored _newWordsComparedTogether: " << _newWordsComparedTogether << std::endl;

	infile->read(reinterpret_cast<char *> (&_lastWordId),sizeof(int)); //_lastWordId
	// std::cout << "restored _lastWordId: " << _lastWordId << std::endl;

	infile->read(reinterpret_cast<char *> (&useDistanceL1_),sizeof(bool)); //useDistanceL1_
	// std::cout << "restored useDistanceL1_: " << useDistanceL1_ << std::endl;

	int strat;
	infile->read(reinterpret_cast<char *> (&strat),sizeof(int)); //useDistanceL1_
	_strategy = static_cast<NNStrategy>(strat);
	// std::cout << "restored _strategy: " << _strategy << std::endl;

	//std::map<int, int> _mapIndexId;
	_mapIndexId.clear();
	int tempa, tempb;
	int mapIndexId_size;
	infile->read(reinterpret_cast<char *>(&mapIndexId_size), sizeof(int));
	int* mem_ptr = (int*) malloc(mapIndexId_size * sizeof(int) * 2);
	int* mem_ptr_temp = mem_ptr;
	int* mem_ptr_temp2 = mem_ptr + 1;
	infile->read(reinterpret_cast<char *>(mem_ptr), mapIndexId_size * sizeof(int) * 2);
	for(int i=0; i<mapIndexId_size; i++)
	{
		_mapIndexId.insert(std::pair<int,int>(*mem_ptr_temp, *mem_ptr_temp2));
		mem_ptr_temp += 2;
		mem_ptr_temp2 += 2;
	}
	// free(mem_ptr);

	// for(int i=0; i<mapIndexId_size; i++)
	// {
	// 	infile->read(reinterpret_cast<char*> (&tempa),sizeof(int));
	// 	infile->read(reinterpret_cast<char*> (&tempb),sizeof(int));
	// 	_mapIndexId.insert(std::pair<int,int>(tempa,tempb));
	// }

	// std::map<int, int> _mapIdIndex;
	_mapIdIndex.clear();
	int mapIdIndex_size;
	infile->read(reinterpret_cast<char *>(&mapIdIndex_size), sizeof(int));
	// mem_ptr = (int*) malloc(mapIdIndex_size * sizeof(int) * 2);
	mem_ptr_temp = mem_ptr;
	mem_ptr_temp2 = mem_ptr + 1;
	infile->read(reinterpret_cast<char *>(mem_ptr), mapIdIndex_size * sizeof(int) * 2);
	// std::cout << "mapidindex size: " << mapIdIndex_size << std::endl;
	// std::cout << "mapIndexId size: " << mapIndexId_size << std::endl;

	for(int i=0; i<mapIdIndex_size; i++)
	{
		_mapIdIndex.insert(std::pair<int,int>(*mem_ptr_temp, *mem_ptr_temp2));
		mem_ptr_temp += 2;
		mem_ptr_temp2 += 2;
	}
	free(mem_ptr);

	// for(int i=0; i<mapIdIndex_size; i++)
	// {
	// 	infile->read(reinterpret_cast<char*> (&tempa),sizeof(int));
	// 	infile->read(reinterpret_cast<char*> (&tempb),sizeof(int));
	// 	_mapIdIndex.insert(std::pair<int,int>(tempa,tempb));
	// }

	// std::map<int, VisualWord *> _unusedWords; //<id,VisualWord*>, note that these words stay in _visualWords
	_unusedWords.clear();
	int unusedwords_size;
	infile->read(reinterpret_cast<char *>(&unusedwords_size), sizeof(int));  //read in the number of visual words
	// std::cout << "restored unused word size: " << unusedwords_size << std::endl;
	for(int i=0; i<unusedwords_size; i++) 
	{
		cv::Mat d;	
		d = cv::Mat(1, 32, CV_8U);	
		infile->read(reinterpret_cast<char *>(&id), sizeof(int));  //read in the id
		v = new VisualWord(id,d);
		v->load_visualword(infile);
		_unusedWords.insert(std::pair<int,VisualWord*>(id,v));
	}

	// std::set<int> _notIndexedWords;			  // Words that are not indexed in the dictionary
	_notIndexedWords.clear();
	int notIndexedWords_size;
	infile->read(reinterpret_cast<char *>(&notIndexedWords_size), sizeof(int));
	// std::cout << "restored notIndexed size: " << notIndexedWords_size << std::endl;
	for(int i=0; i<notIndexedWords_size; i++)
	{
		infile->read(reinterpret_cast<char*> (&tempa),sizeof(int));
		_notIndexedWords.insert(tempa);
	}

	// std::set<int> _removedIndexedWords;		  // Words not anymore in the dictionary but still indexed in the dictionary
	_removedIndexedWords.clear();
	int removedIndexedWords_size;
	infile->read(reinterpret_cast<char *>(&removedIndexedWords_size), sizeof(int));
	// std::cout << "restored removeIndexed size: " << removedIndexedWords_size << std::endl;
	for(int i=0; i<removedIndexedWords_size; i++)
	{
		infile->read(reinterpret_cast<char*> (&tempa),sizeof(int));
		_removedIndexedWords.insert(tempa);
	}

	infile->close();

	_dataTree = cv::Mat();

	auto t3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms2 = t3 - t2;
	std::cout << "vwdictionary load time 2 in milliseconds: " << fp_ms2.count() << std::endl;

	//load the flann index data
	_flannIndex = new FlannIndex();
	_flannIndex->load_index();
}

} // namespace rtabmap
