#pragma once


#include <exception>
#include <string>
#include <utility>

class FileNotFoundException : public std::exception {
private:
    std::string file_path;
public:
    FileNotFoundException(std::string path) : file_path(std::move(path)) {}

    std::string what() {
        return file_path;
    }
};
