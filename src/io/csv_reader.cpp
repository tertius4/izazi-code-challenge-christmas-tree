#include "io/csv_reader.h"

#include <cstddef>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace packer::io {

std::vector<std::vector<std::string>> readCsv(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open " + path);

    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    std::vector<std::vector<std::string>> rows;
    std::vector<std::string> row;
    std::string field;
    bool inQuotes = false;
    bool fieldStarted = false; // true once the current row has any content

    size_t i = 0;
    if (text.compare(0, 3, "\xEF\xBB\xBF") == 0) i = 3; // skip UTF-8 BOM

    auto endField = [&] {
        row.push_back(std::move(field));
        field.clear();
    };

    auto endRow = [&] {
        if (fieldStarted) { // skip blank lines
            endField();
            rows.push_back(std::move(row));
            row.clear();
        }
        fieldStarted = false;
    };

    for (; i < text.size(); ++i) {
        const char c = text[i];

        if (inQuotes) {
            if (c == '"') {
                const bool doubledQuote = i + 1 < text.size() && text[i + 1] == '"';
                if (doubledQuote) {
                    field += '"'; // stands for one literal quote
                    ++i;
                } else {
                    inQuotes = false;
                }
            } else {
                field += c;
            }
            continue;
        }

        switch (c) {
        case '"':
            inQuotes = true;
            fieldStarted = true;
            break;
        case ',':
            endField();
            fieldStarted = true;
            break;
        case '\r':
            break; // ignore; the following '\n' ends the row
        case '\n':
            endRow();
            break;
        default:
            field += c;
            fieldStarted = true;
        }
    }

    if (inQuotes) throw std::runtime_error("Unclosed quote in " + path);
    endRow(); // last line without a trailing newline

    return rows;
}

} // namespace packer::io
