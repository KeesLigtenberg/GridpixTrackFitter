#include "makeNoisyPixelMask.h"

#include <iostream>

#include "TH2.h"
//#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "getHistFromTree.h"

//threshold is maximum of number times mean
pixelMask makeNoisyPixelMask(TTree* hitTable, int plane, double threshold, std::pair<int,int> gridsize ) {
	int gridx=gridsize.first, gridy=gridsize.second;
	auto p=std::to_string(plane);
	TH2D hist("noisyPixelMaskHistogram", "histogram of all pixels",gridx,0,gridx,gridy,0,gridy );
	getHistFromTree(*hitTable, "mimosa["+p+"][].row:mimosa["+p+"][].column", "1", hist.GetName(),"goff"/*, DEBUG 1e5 */);

	auto mean = hist.GetEntries()/gridx/gridy;
//	std::cout<<"mean is "<<mean<<std::endl;
	pixelMask mask(gridx, std::vector<char>(gridy, 0) );
	int nmasked=0;
	for(int x=0; x<gridx; x++) {
		for(int y=0; y<gridy; y++) {
			if( hist.GetBinContent(x+1, y+1) > mean*threshold ) {
//				std::cout<<hist.GetBinContent(x+1, y+1)<<" > "<<mean*threshold<<std::endl;
				mask[x][y]=1;
				++nmasked;
			}
		}
	}

	std::cout<<"masked "<<nmasked<<" pixels for plane "<<plane<<std::endl;

	return mask;
}

std::vector<Hit> applyPixelMask(const pixelMask& mask, const std::vector<Hit>& hv ) {
	std::vector<Hit> newhv;
	for(const Hit& h:hv)	{
		if(!mask.at(h.column).at(h.row)) {
			newhv.emplace_back(h);
		}
	}
//	std::cout<<"returning hv with size " << newhv.size()<<std::endl;
	return newhv;
}

