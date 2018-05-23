//rootmacro
#ifndef GETOBJECTSFROMFILE_H
#define GETOBJECTSFROMFILE_H
#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

template<class T>
T* getObjectFromFile(std::string objectname, TFile* file) {
	T* obj=dynamic_cast<T*>(file->Get(objectname.c_str()));
	if(!obj) {
		std::cerr<<"could not get object "<<objectname<<" from file"<<std::endl;
		throw 1;	
	}
	return obj;
}

TFile* openFile(std::string filename){ 
	TFile* file = TFile::Open(filename.c_str());
	if(!file) {
		std::cerr<<"could not open file "<<filename<<std::endl;
		throw 1;
	}
	return file;
}

template<class T>
T* getObjectFromFile(std::string objectname, std::string filename) {
	TFile* file=openFile(filename);
	return getObjectFromFile<T>(objectname, file);
}

#endif
