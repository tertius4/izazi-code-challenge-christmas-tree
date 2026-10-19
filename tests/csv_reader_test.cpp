#include "io/csv_reader.h"

#include <stdexcept>

#include "temp_file.h"

using packer::io::readCsv;
using Rows = std::vector<std::vector<std::string>>;

class CsvReaderTest : public TempDirTest {};

TEST_F(CsvReaderTest, ParsesSimpleRows) {
    const auto file = writeFile("a.csv", "a,b,c\n1,2,3\n");
    EXPECT_EQ(readCsv(file), (Rows{{"a", "b", "c"}, {"1", "2", "3"}}));
}

TEST_F(CsvReaderTest, QuotedFieldKeepsComma) {
    const auto file = writeFile("a.csv", "\"x,y\",z\n");
    EXPECT_EQ(readCsv(file), (Rows{{"x,y", "z"}}));
}

TEST_F(CsvReaderTest, DoubledQuotesBecomeOneQuote) {
    const auto file = writeFile("a.csv", "\"say \"\"hi\"\"\",b\n");
    EXPECT_EQ(readCsv(file), (Rows{{"say \"hi\"", "b"}}));
}

TEST_F(CsvReaderTest, NewlineInsideQuotesStaysInField) {
    const auto file = writeFile("a.csv", "\"line1\nline2\",b\nc,d\n");
    EXPECT_EQ(readCsv(file), (Rows{{"line1\nline2", "b"}, {"c", "d"}}));
}

TEST_F(CsvReaderTest, CrlfLineEndings) {
    const auto file = writeFile("a.csv", "a,b\r\n1,2\r\n");
    EXPECT_EQ(readCsv(file), (Rows{{"a", "b"}, {"1", "2"}}));
}

TEST_F(CsvReaderTest, StripsUtf8Bom) {
    const auto file = writeFile("a.csv", "\xEF\xBB\xBF"
                                         "a,b\n");
    EXPECT_EQ(readCsv(file), (Rows{{"a", "b"}}));
}

TEST_F(CsvReaderTest, SkipsBlankLines) {
    const auto file = writeFile("a.csv", "a\n\n\nb\n");
    EXPECT_EQ(readCsv(file), (Rows{{"a"}, {"b"}}));
}

TEST_F(CsvReaderTest, PreservesRowsOfDifferentLengths) {
    const auto file = writeFile("a.csv", "a,b,c\n1\n1,2\n");
    EXPECT_EQ(readCsv(file), (Rows{{"a", "b", "c"}, {"1"}, {"1", "2"}}));
}

TEST_F(CsvReaderTest, LastRowWithoutTrailingNewline) {
    const auto file = writeFile("a.csv", "a,b\n1,2");
    EXPECT_EQ(readCsv(file), (Rows{{"a", "b"}, {"1", "2"}}));
}

TEST_F(CsvReaderTest, EmptyFieldsArePreserved) {
    const auto file = writeFile("a.csv", "a,,c\n,,\n");
    EXPECT_EQ(readCsv(file), (Rows{{"a", "", "c"}, {"", "", ""}}));
}

TEST_F(CsvReaderTest, EmptyFileHasNoRows) {
    const auto file = writeFile("a.csv", "");
    EXPECT_TRUE(readCsv(file).empty());
}

TEST_F(CsvReaderTest, MissingFileThrows) {
    EXPECT_THROW(readCsv(path("does_not_exist.csv")), std::runtime_error);
}

TEST_F(CsvReaderTest, UnclosedQuoteThrows) {
    const auto file = writeFile("a.csv", "\"never closed,b\n");
    EXPECT_THROW(readCsv(file), std::runtime_error);
}

TEST(CsvReaderRealData, ParsesTreeCsv) {
    const auto rows = readCsv(TREE_CSV_PATH);
    ASSERT_EQ(rows.size(), 16u); // header + 15 vertices
    EXPECT_EQ(rows[0], (std::vector<std::string>{"vertex", "x", "y", "description"}));
    for (const auto &row : rows) {
        EXPECT_EQ(row.size(), 4u);
    }
}
