include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(tinyxml)

set(SRC
	${tinyxml_SRC}
	${CMAKE_CURRENT_SOURCE_DIR}/cyBVH.h
	${CMAKE_CURRENT_SOURCE_DIR}/cyColor.h
	${CMAKE_CURRENT_SOURCE_DIR}/cyCore.h
	${CMAKE_CURRENT_SOURCE_DIR}/cyMatrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/cyPoint.h
	${CMAKE_CURRENT_SOURCE_DIR}/cyTriMesh.h
	${CMAKE_CURRENT_SOURCE_DIR}/lights.h
	${CMAKE_CURRENT_SOURCE_DIR}/lodepng.h
	${CMAKE_CURRENT_SOURCE_DIR}/lodepng.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/materials.h
	${CMAKE_CURRENT_SOURCE_DIR}/objects.h
	${CMAKE_CURRENT_SOURCE_DIR}/scene.h
	${CMAKE_CURRENT_SOURCE_DIR}/texture.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/texture.h
	${CMAKE_CURRENT_SOURCE_DIR}/viewport.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/xmlload.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
	)

set(project8_data_gathering_SRC 
	${SRC}
	PARENT_SCOPE)