# GridpixTrackFitter
Track fitting and analysis code for the data taken with a single chip gridpix detector in July 2017

To be run with Root cern version 6.08.02 and at least C++11

Use FitTracksTimePix.cpp to fit tracks to timepix: 

	root
	[] .X FitTracksTimePix.cpp+("example_converted.root") 

In root the compilation '+' is required. Note that some functions are in unusual places, such as the many static functions in HoughTransformer. To code can be cleaned up a lot.

To convert the raw tree to a more convenient format, I used my convert script. There is also some code there to find the trigger offset, needed for synchronising it with the telescope.

The description here is very limited: 

please contact me at: *cligtenb(at)nikhef.nl*
