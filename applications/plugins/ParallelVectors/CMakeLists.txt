cmake_minimum_required(VERSION 3.1)
project(ParallelVectors)

set(HEADER_FILES
    initParallelVectors.h
    vector_parallel.h
    ParallelMemoryManager.h
    ParallelTypes.h
#    ParallelState.h
#    ParallelMechanicalState.h
    ParallelMechanicalObject.h
)

set(SOURCE_FILES
    initParallelVectors.cpp
    ParallelTypes.cpp
#    ParallelState.cpp
#    ParallelMechanicalState.cpp
    ParallelMechanicalObject.cpp
    ParallelMechanicalObject.inl
)

find_package(SofaFramework REQUIRED)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_PARALLELVECTORS")
target_link_libraries(${PROJECT_NAME} SofaCore SofaBaseMechanics)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/..")

install(TARGETS ${PROJECT_NAME}
        COMPONENT ParallelVectors_libraries
        EXPORT ParallelVectorsTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
