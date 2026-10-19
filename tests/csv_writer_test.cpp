#include "io/csv_writer.h"

#include "io/csv_reader.h"
#include "temp_file.h"

using packer::io::readCsv;
using packer::io::writeCsv;
using Rows = std::vector<std::vector<std::string>>;

class CsvWriterTest : public TempDirTest {};

TEST_F(CsvWriterTest, WritesPlainFields) {
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, Rows{{"a", "b"}, {"1", "2"}}));
    EXPECT_EQ(slurp(file), "a,b\n1,2\n");
}

TEST_F(CsvWriterTest, QuotesFieldWithComma) {
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, Rows{{"x,y", "z"}}));
    EXPECT_EQ(slurp(file), "\"x,y\",z\n");
}

TEST_F(CsvWriterTest, EscapesQuotesByDoubling) {
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, Rows{{"say \"hi\""}}));
    EXPECT_EQ(slurp(file), "\"say \"\"hi\"\"\"\n");
}

TEST_F(CsvWriterTest, QuotesFieldWithNewline) {
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, Rows{{"a\nb", "c"}}));
    EXPECT_EQ(slurp(file), "\"a\nb\",c\n");
}

TEST_F(CsvWriterTest, WritesEmptyFields) {
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, Rows{{"a", "", "c"}}));
    EXPECT_EQ(slurp(file), "a,,c\n");
}

TEST_F(CsvWriterTest, EmptyDataWritesEmptyFile) {
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, Rows{}));
    EXPECT_EQ(slurp(file), "");
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
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, original));
    EXPECT_EQ(readCsv(file), original);
}

TEST_F(CsvWriterTest, RoundTripsRealData) {
    const auto original = readCsv(TREE_CSV_PATH);
    const auto file = path("out.csv");
    ASSERT_TRUE(writeCsv(file, original));
    EXPECT_EQ(readCsv(file), original);
}
