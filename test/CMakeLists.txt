find_package (Catch2 REQUIRED)
find_package(rdkit REQUIRED)


add_executable(unittests
        test_runner.cpp
        )

include_directories(.)

target_link_libraries(unittests PRIVATE
        ${PROJECT_NAME}
        ${CONAN_LIBS}
        )

add_test(NAME Test COMMAND unittests)
