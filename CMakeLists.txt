PROJECT(Nanopositioning_SEM_win)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#FIND_PACKAGE( OpenCV REQUIRED )

FIND_PACKAGE(VISP REQUIRED)
IF(VISP_FOUND)
  INCLUDE(${VISP_USE_FILE})
ENDIF(VISP_FOUND)


#INCLUDE_DIRECTORIES(test)

ADD_EXECUTABLE( testImageAcquisition 
		testImageAcquisition.cpp
		) 
ADD_EXECUTABLE( testRobotControl 
		testRobotControl.cpp
		) 


ADD_EXECUTABLE( Nanopositioning_SEM_win 
		Nanopositioning_SEM_win.cpp
 		npFeatureLuminance.cpp 
		npFeatureLuminance.h  
		npRegression.cpp 
		npRegression.h
		npImageProcessing.cpp
		npImageProcessing.h
		) 
