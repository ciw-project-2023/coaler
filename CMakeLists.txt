cmake_minimum_required(VERSION 3.15..3.27)

project(aligner CXX)

find_package(rdkit REQUIRED)

set(SOURCE_FILES
  src/main.cpp
  src/SingleAligner.cpp)


# Outputs compile commands used by clang-tidy.
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})

add_executable(${PROJECT_NAME} STATIC ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})

target_include_directories(${PROJECT_NAME} PUBLIC
            ${CMAKE_BINARY_DIR}/generated/
            )

target_compile_options(${PROJECT_NAME} PRIVATE
            -Wall
            -Wextra
            )

# Testing
enable_testing()
add_subdirectory (test)
