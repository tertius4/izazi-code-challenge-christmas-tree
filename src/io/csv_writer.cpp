#include "csv_writer.h"

#include <fstream>

// Wraps a field in quotes if it contains a comma, quote or newline.
static std::string escapeField(const std::string &field)
{
    if (field.find_first_of(",\"\n") == std::string::npos)
    {
        return field;
    }
    std::string result = "\"";
    for (char c : field)
    {
        if (c == '"')
            result += "\"\""; // a quote is escaped by doubling it
        else
            result += c;
    }
    result += "\"";
    return result;
}

namespace packer::io
{
    bool writeCsv(const std::string &path, const std::vector<std::vector<std::string>> &data)
    {
        std::ofstream out(path);
        if (!out)
            return false;

        for (const auto &row : data)
        {
            for (size_t i = 0; i < row.size(); ++i)
            {
                if (i > 0)
                    out << ",";
                out << escapeField(row[i]);
            }
            out << "\n";
        }
        return out.good();
    }
}