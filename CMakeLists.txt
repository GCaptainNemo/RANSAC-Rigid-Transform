cmake_minimum_required(VERSION 3.0)

project(RANSAC_RIGID_TRANSFORM)


find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(include ${OpenCV_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})


add_executable(RANSAC_RIGID_TRANSFORM src/main.cpp)


target_link_libraries(RANSAC_RIGID_TRANSFORM ${OpenCV_LIBRARIES})

