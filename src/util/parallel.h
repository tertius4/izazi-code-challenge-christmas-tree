#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace packer::util {
// Every core of the machine, or 1 if that is unknown.
inline int allCores() {
    return static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
}

// Calls body(i) for every i in [0, count) on `threads` threads (0: every core). Each thread
// takes the next unclaimed i as soon as it is free, so slow and quick iterations even out.
// Returns once every call has finished; if any threw, rethrows the first exception.
template <class Body> void parallelFor(std::size_t count, int threads, Body body) {
    if (count == 0) return;
    const std::size_t want = static_cast<std::size_t>(threads > 0 ? threads : allCores());
    const std::size_t workers = std::min(want, count);

    std::atomic<std::size_t> next{0};
    std::exception_ptr error;
    std::mutex errorMutex;
    auto work = [&] {
        for (std::size_t i = next++; i < count; i = next++) {
            try {
                body(i);
            } catch (...) {
                const std::lock_guard lock(errorMutex);
                if (!error) error = std::current_exception();
                next = count; // the others finish their current call and stop
            }
        }
    };

    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < workers; ++t) {
        pool.emplace_back(work);
    }
    work(); // the calling thread works too
    for (std::thread &t : pool) {
        t.join();
    }

    if (error) std::rethrow_exception(error);
}
} // namespace packer::util
