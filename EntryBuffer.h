/*
 * EntryBuffer.h
 *
 *  Created on: Oct 4, 2017
 *      Author: cligtenb
 *
 *      Specialised for use with tpcEntryNumbers
 *
 *      For other purposes the container (now map) should be reconsidered.
 */

#ifndef ENTRYBUFFER_H_
#define ENTRYBUFFER_H_
#include <map>
#include <functional>

template <class K, class V>
class EntryBuffer {
public:
	void placeInBuffer(K tpcEntryNumber, const V&);
	bool isInBuffer(K tpcEntryNumber);
	V& getFromBuffer(K tpcEntryNumber);
	void removeFromBuffer(K tpcEntryNumber);
	void writeBufferUpTo(K tpcEntryNumber, std::function<void(V&)> function);
	void writeBuffer(std::function<void(V&)> function);
private:
	std::map<K,V> buffer; //buffer by tpcEntryNumber
};


template<class K, class V>
inline void EntryBuffer<K, V>::placeInBuffer(K tpcEntryNumber, const V& entry) {
	if(not buffer.insert( {tpcEntryNumber, entry} ).second) {
		throw "Unexpected element in buffer!";
	}
}

template<class K, class V>
inline bool EntryBuffer<K, V>::isInBuffer(K tpcEntryNumber) {
	if( buffer.find(tpcEntryNumber) == buffer.end() ) {
		return false;
	} else {
		return true;
	}
}

template<class K, class V>
inline V& EntryBuffer<K, V>::getFromBuffer(K tpcEntryNumber) {
	return buffer.at(tpcEntryNumber);
}

template <class K, class V>
inline void EntryBuffer<K,V>::removeFromBuffer(K tpcEntryNumber) {
	auto result=buffer.find(tpcEntryNumber);
	if(result==buffer.end()) throw "Could not find entry in buffer";
	buffer.erase(result);
}

template <class K, class V>
inline void EntryBuffer<K,V>::writeBufferUpTo(K tpcEntryNumber, std::function<void(V&)> writeFunction) {
	for(auto e=buffer.begin();
			 e!=buffer.end() and e->first < tpcEntryNumber;
			 e=buffer.erase(e) ) {
		writeFunction(e->second);
	}
}

template <class K, class V>
void EntryBuffer<K,V>::writeBuffer(std::function<void(V&)> writeFunction) {
	for(auto& e : buffer) {
		writeFunction(e.second);
	}
	buffer.clear();
}



#endif /* ENTRYBUFFER_H_ */
