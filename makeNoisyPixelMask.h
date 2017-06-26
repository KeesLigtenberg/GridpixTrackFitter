#ifndef MAKENOISYPIXELMASK_H
#define MAKENOISYPIXELMASK_H

#include <vector>
#include <utility>
#include <string>

#include "TTree.h"

#include "Hit.h"

using pixelMask=std::vector< std::vector<char> >;

//threshold is maximum of number times mean
pixelMask makeNoisyPixelMask(TTree* hitTable, int plane, double threshold, std::pair<int,int> gridsize={1152,576} ) ;

std::vector<Hit> applyPixelMask(const pixelMask& mask, const std::vector<Hit>& hv ) ;

#endif
