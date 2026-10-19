#include "io/csv_writer.h"

#include <cstddef>
#include <fstream>

// Wraps a field in quotes if it contains a comma, quote or line break.
static std::string escapeField(const std::string &field) {
    if (field.find_first_of(",\"\r\n") == std::string::npos) {
        return field;
    }

    std::string result = "\"";
    for (const char ch : field) {
        if (ch == '"') {
            result += "\"\""; // a quote is escaped by doubling it
        } else {
            result += ch;
        }
    }

    result += "\"";

    return result;
}

namespace packer::io {
bool writeCsv(const std::string &path, const std::vector<std::vector<std::string>> &rows) {
    std::ofstream out(path);
    if (!out) return false;

    for (const auto &row : rows) {
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i > 0) out << ",";
            out << escapeField(row[i]);
        }

        out << "\n";
    }

    return out.good();
}
} // namespace packer::io
