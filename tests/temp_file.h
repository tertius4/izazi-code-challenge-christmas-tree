#pragma once

#include <atomic>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include <gtest/gtest.h>

// Gives each test its own scratch directory, removed afterwards.
class TempDirTest : public ::testing::Test {
  protected:
    void SetUp() override {
        static std::atomic<int> counter{0};
        const auto *info = ::testing::UnitTest::GetInstance()->current_test_info();
        dir_ = std::filesystem::temp_directory_path() /
               ("packer_tests_" + std::string(info->test_suite_name()) + "_" +
                std::string(info->name()) + "_" + std::to_string(counter++));
        std::filesystem::create_directories(dir_);
    }

    void TearDown() override {
        std::error_code ec;
        std::filesystem::remove_all(dir_, ec);
    }

    std::string path(const std::string &name) const {
        return (dir_ / name).string();
    }

    std::string writeFile(const std::string &name, const std::string &content) const {
        const std::string p = path(name);
        std::ofstream out(p, std::ios::binary);
        out << content;
        return p;
    }

    static std::string slurp(const std::string &p) {
        std::ifstream in(p, std::ios::binary);
        return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    }

    std::filesystem::path dir_;
};
