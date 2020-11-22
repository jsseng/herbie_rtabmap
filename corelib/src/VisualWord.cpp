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

#include "rtabmap/core/VisualWord.h"
#include "rtabmap/utilite/ULogger.h"
#include "rtabmap/utilite/UStl.h"

using std::cout;

namespace rtabmap
{

VisualWord::VisualWord(int id, const cv::Mat & descriptor, int signatureId) :
	_id(id),
	_descriptor(descriptor),
	_saved(false),
	_totalReferences(0)
{
	if(signatureId)
	{
		addRef(signatureId);
	}
}

VisualWord::~VisualWord()
{
}

void VisualWord::addRef(int signatureId)
{
	std::map<int, int>::iterator iter = _references.find(signatureId);
	if(iter != _references.end())
	{
		(*iter).second += 1;
	}
	else
	{
		_references.insert(std::pair<int, int>(signatureId, 1));
	}
	++_totalReferences;
}

int VisualWord::removeAllRef(int signatureId)
{
	int removed = uTake(_references, signatureId, 0);
	_totalReferences -= removed;
	return removed;
}

unsigned long VisualWord::getMemoryUsed() const
{
	unsigned long memoryUsage = sizeof(VisualWord);
	memoryUsage += _references.size() * (sizeof(int)*2+sizeof(std::map<int ,int>::iterator)) + sizeof(std::map<int ,int>);
	memoryUsage += _oldReferences.size() * (sizeof(int)*2+sizeof(std::map<int ,int>::iterator)) + sizeof(std::map<int ,int>);
	memoryUsage += _descriptor.total() * _descriptor.elemSize();
	return memoryUsage;
}

void VisualWord::debug() {
	std::cout << "_references: " << _references.size() << "," << _oldReferences.size() << std::endl;
}

void VisualWord::save_visualword(std::ofstream *outfile) {
	// int data = _descriptor.cols;
	int a,b;
	outfile->write(reinterpret_cast<const char *> (_descriptor.data),32); //32 bytes = 256 ORB bits / 8 bits

	outfile->write(reinterpret_cast<const char *> (&_saved),sizeof(bool)); //_saved

	outfile->write(reinterpret_cast<const char *> (&_totalReferences),sizeof(int)); //_totalReferences

	int _references_size = _references.size();
	outfile->write(reinterpret_cast<const char *> (&_references_size),sizeof(int)); //_references size
	for(std::map<int, int>::iterator iter = _references.begin(); iter!=_references.end(); ++iter)
	{
		a = iter->first;
		b = iter->second;
		outfile->write(reinterpret_cast<char*> (&a),sizeof(int));
		outfile->write(reinterpret_cast<char*> (&b),sizeof(int));
	}

	int _oldReferences_size = _oldReferences.size();
	outfile->write(reinterpret_cast<const char *> (&_oldReferences_size),sizeof(int)); //_references size
	for(std::map<int, int>::iterator iter = _oldReferences.begin(); iter!=_oldReferences.end(); ++iter)
	{
		a = iter->first;
		b = iter->second;
		outfile->write(reinterpret_cast<char*> (&a),sizeof(int));
		outfile->write(reinterpret_cast<char*> (&b),sizeof(int));
	}
}

void VisualWord::load_visualword(std::ifstream *infile) {
	// int data = _descriptor.cols;
	int a,b;
	infile->read(reinterpret_cast<char *> (_descriptor.data),32); //32 bytes = 256 ORB bits / 8 bits

	infile->read(reinterpret_cast<char *> (&_saved),sizeof(bool)); //_saved

	infile->read(reinterpret_cast<char *> (&_totalReferences),sizeof(int)); //_totalReferences

	int _references_size;
	infile->read(reinterpret_cast<char *> (&_references_size),sizeof(int)); //_references size
	for (int i=0; i<_references_size; i++)
	{
		infile->read(reinterpret_cast<char*> (&a),sizeof(int));
		infile->read(reinterpret_cast<char*> (&b),sizeof(int));
		_references.insert(std::pair<int,int>(a,b));
	}

	int _oldReferences_size;
	infile->read(reinterpret_cast<char *> (&_oldReferences_size),sizeof(int)); //_oldReferences size
	for (int i=0; i<_oldReferences_size; i++)
	{
		infile->read(reinterpret_cast<char*> (&a),sizeof(int));
		infile->read(reinterpret_cast<char*> (&b),sizeof(int));
		_oldReferences.insert(std::pair<int,int>(a,b));
	}
}

} // namespace rtabmap
