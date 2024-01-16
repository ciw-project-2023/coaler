#pragma once

#include <exception>
#include <string>
#include <utility>

namespace coaler::io {
    class FileNotFoundException : public std::exception {
      private:
        std::string m_file_path;

      public:
        explicit FileNotFoundException(std::string path) : m_file_path(std::move(path)) {}

        std::string what() { return m_file_path; }
    };
}  // namespace coaler::io
