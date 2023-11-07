//
// Created by niklas on 11/7/23.
//

#ifndef ALIGNER_FILENOTFOUNDEXCEPTION_H
#define ALIGNER_FILENOTFOUNDEXCEPTION_H


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


#endif //ALIGNER_FILENOTFOUNDEXCEPTION_H
