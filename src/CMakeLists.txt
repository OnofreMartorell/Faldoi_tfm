project(faldoi_qtcreator)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)


set(CMAKE_BUILD_TYPE DEBUG)

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})

find_package(JPEG REQUIRED)
include_directories(${JPEG_INCLUDE_DIR})

find_package(TIFF REQUIRED)
include_directories(${TIFF_INCLUDE_DIR})

find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3 -std=c99")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -O0 -ggdb -DNDEBUG -Wall -Wextra -std=c99")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -ggdb -DNDEBUG -Wall -Wextra -std=c++11")

SET(SHARED_C_SRC iio.c mask.c xmalloc.c bicubic_interpolation.c elap_recsep.c)
SET(SHARED_CPP_SRC tvl2_model.cpp nltv_model.cpp tvcsad_model.cpp nltvcsad_model.cpp tvl2w_model.cpp nltvcsadw_model.cpp nltvw_model.cpp tvcsadw_model.cpp aux_energy_model.cpp energy_model.cpp heuristic_interpolation.cpp)

message(${PNG_LIBRARY})
add_executable(sparse_flow ${SHARED_C_SRC} ${SHARED_CPP_SRC} sparse_flow.cpp)
add_executable(local_faldoi ${SHARED_C_SRC} ${SHARED_CPP_SRC} local_faldoi.cpp)
add_executable(global_faldoi ${SHARED_C_SRC} ${SHARED_CPP_SRC} global_faldoi.cpp)
target_link_libraries(sparse_flow ${PNG_LIBRARY} ${JPEG_LIBRARY} ${TIFF_LIBRARY})
target_link_libraries(local_faldoi ${PNG_LIBRARY} ${JPEG_LIBRARY} ${TIFF_LIBRARY})
target_link_libraries(global_faldoi ${PNG_LIBRARY} ${JPEG_LIBRARY} ${TIFF_LIBRARY})
