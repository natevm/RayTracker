include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set (tinyxml_SRC 
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
	${CMAKE_CURRENT_SOURCE_DIR}/tinystr.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/tinystr.h
	${CMAKE_CURRENT_SOURCE_DIR}/tinyxml.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tinyxml.h
	${CMAKE_CURRENT_SOURCE_DIR}/tinyxmlerror.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tinyxmlparser.cpp
	PARENT_SCOPE)
