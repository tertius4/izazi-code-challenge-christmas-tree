#pragma once

#include <string>
#include <vector>

namespace packer::io
{
    // Writes rows of strings to a CSV file. Returns true on success.
    bool writeCsv(const std::string &path, const std::vector<std::vector<std::string>> &data);
}