#pragma once

#include <exception>
#include <string>
#include <utility>

namespace coaler {
    namespace io {
        class FileNotFoundException : public std::exception {
          private:
            std::string m_file_path;

          public:
            FileNotFoundException(std::string path) : m_file_path(std::move(path)) {}

            const char* what() const noexcept override { return m_file_path.c_str(); }
        };
    }  // namespace io
}  // namespace coaler
