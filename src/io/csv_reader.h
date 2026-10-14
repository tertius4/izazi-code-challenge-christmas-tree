#pragma once

#include <string>
#include <vector>

namespace packer::io
{
    // Reads a CSV file into rows of string fields.
    // Handles quoted fields (commas, doubled quotes and newlines inside quotes),
    // both LF and CRLF line endings, and a leading UTF-8 BOM.
    // Empty lines are skipped. Rows may have different lengths.
    // Throws std::runtime_error if the file cannot be opened or a quote is unclosed.
    std::vector<std::vector<std::string>> readCsv(const std::string &path);

}