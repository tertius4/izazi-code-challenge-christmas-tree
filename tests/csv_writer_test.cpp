#include "io/csv_writer.h"

#include "io/csv_reader.h"
#include "temp_file.h"

using packer::io::readCsv;
using packer::io::writeCsv;
using Rows = std::vector<std::vector<std::string>>;

class CsvWriterTest : public TempDirTest {};

TEST_F(CsvWriterTest, WritesPlainFields) {
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, Rows{{"a", "b"}, {"1", "2"}}));
    EXPECT_EQ(slurp(p), "a,b\n1,2\n");
}

TEST_F(CsvWriterTest, QuotesFieldWithComma) {
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, Rows{{"x,y", "z"}}));
    EXPECT_EQ(slurp(p), "\"x,y\",z\n");
}

TEST_F(CsvWriterTest, EscapesQuotesByDoubling) {
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, Rows{{"say \"hi\""}}));
    EXPECT_EQ(slurp(p), "\"say \"\"hi\"\"\"\n");
}

TEST_F(CsvWriterTest, QuotesFieldWithNewline) {
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, Rows{{"a\nb", "c"}}));
    EXPECT_EQ(slurp(p), "\"a\nb\",c\n");
}

TEST_F(CsvWriterTest, WritesEmptyFields) {
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, Rows{{"a", "", "c"}}));
    EXPECT_EQ(slurp(p), "a,,c\n");
}

TEST_F(CsvWriterTest, EmptyDataWritesEmptyFile) {
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, Rows{}));
    EXPECT_EQ(slurp(p), "");
}

TEST_F(CsvWriterTest, ReturnsFalseWhenPathIsUnwritable) {
    EXPECT_FALSE(writeCsv(path("no_such_dir/out.csv"), Rows{{"a"}}));
}

TEST_F(CsvWriterTest, RoundTripsAwkwardFields) {
    const Rows original{
        {"plain", "with,comma", "with \"quote\""},
        {"multi\nline", "carriage\rreturn", ""},
        {"crlf\r\nfield", ",", "\""},
    };
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, original));
    EXPECT_EQ(readCsv(p), original);
}

TEST_F(CsvWriterTest, RoundTripsRealData) {
    const auto original = readCsv(TREE_CSV_PATH);
    const auto p = path("out.csv");
    ASSERT_TRUE(writeCsv(p, original));
    EXPECT_EQ(readCsv(p), original);
}
