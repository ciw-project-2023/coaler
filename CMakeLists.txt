cmake_minimum_required(VERSION 3.15..3.27)

option(USE_CONAN "enables the dependency management via conan for this package" ON)

project(aligner CXX)

if (USE_CONAN)
    if (EXISTS ${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
        set(CONAN_SYSTEM_INCLUDES ON)
        include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
        conan_basic_setup()

        # Enable sanitizer and coverage if conan profile propagates them.
        if (${CONAN_SETTINGS_COMPILER_SANITIZER} MATCHES "AddressUndefinedBehavior")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -fno-sanitize=alignment")
        endif ()

        if (${CONAN_SETTINGS_COMPILER_COVERAGE})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
        endif ()
    else ()
        message(WARNING "Dependencies not found. Run \"conan install\" first.")
    endif ()
endif ()

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
