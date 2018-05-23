# GridpixTrackFitter
Track fitting and analysis code for the data taken with a single chip gridpix detector in July 2017

Run the scripts with Cern-ROOT version 6.08.02 and at least C++11. The data can be found at [zenodo](https://doi.org/10.5281/zenodo.1251806).


Note that some functions are in unusual places, such as the many static functions in HoughTransformer. The code can be cleaned up a lot.

##convert from raw
To convert the raw tree to a more convenient format, I used my convert script. There is also some code there to find the trigger offset, needed for synchronising it with the telescope. I recommend to start from the converted files

##Fit track in Timepix only
Use FitTracksTimePix.cpp to fit tracks to timepix: 

	root
	[] .X FitTracksTimePix.cpp+("example_converted.root") 

In root the compilation '+' is required.

##Fit tracks Timepix with telescope
Run CombineTracks.cpp to get an output file *fitResults.root* with all fits and residuals.

	root
	[] .X CombineTracks.cpp+("58_mimosa_telescope_nikhef_m26_telescope_scan_combined.root","W0015_H04-170713-071936-347_converted.root")

##Make histograms
In order to acquire a file with histograms *histograms.root* from *fitResults.root*, run 

	root
	[] .L /project/lepcol/users/cligtenb/testbeamscripts/trackFitter/resultProcessor.cpp+
	[] resultProcessor p
	[] p.Loop();



The description here is very limited, so please contact me at: *cligtenb(at)nikhef.nl*
