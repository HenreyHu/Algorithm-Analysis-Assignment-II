// run_tests.cpp - runs every SC_TEST that the present module headers register.
// Owner: P1. Modules are picked up automatically, so adding one never edits this file.
//
//   run_tests [--filter PREFIX] [--csv PATH]
//
// Prints PASS/FAIL per test, writes id,name,status,ms,message to the CSV
// (default results/tests.csv; a --filter run writes a CSV only when --csv is given, so a partial
// run never replaces the full list), and ends with "ALL N TESTS PASSED" (exit 0)
// or "k OF N TESTS FAILED" (exit 1).
#define SC_TESTS
#include "core.h"
#if __has_include("energy.h")
#include "energy.h"
#endif
#if __has_include("dp.h")
#include "dp.h"
#endif
#if __has_include("greedy.h")
#include "greedy.h"
#endif
#if __has_include("graphcut.h")
#include "graphcut.h"
#endif
#if __has_include("carve.h")
#include "carve.h"
#endif

#include <optional>
#include <thread>

int main(int argc, char** argv) {
    sc::PauseOnExit pause_on_exit;  // a double-clicked window must not close before the result can be read
    std::string filter, csv_path;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--filter" && i + 1 < argc)
            filter = argv[++i];
        else if (a == "--csv" && i + 1 < argc)
            csv_path = argv[++i];
        else {
            std::fprintf(stderr, "usage: run_tests [--filter PREFIX] [--csv PATH]\n");
            return 2;
        }
    }
    if (csv_path.empty() && filter.empty()) csv_path = "results/tests.csv";

    // Smoke check of the toolchain's thread support, which the energy and DP modules rely on.
    {
        int x = 0;
        std::thread t([&] { x = 1; });
        t.join();
        if (x != 1) {
            std::fprintf(stderr, "std::thread smoke check failed\n");
            return 1;
        }
    }

    std::filesystem::create_directories("tests_tmp");
    std::optional<sc::Csv> csv;
    if (!csv_path.empty()) csv.emplace(csv_path, std::vector<std::string>{"id", "name", "status", "ms", "message"});
    int run = 0, failed = 0;
    for (const sc::TestCase& t : sc::test_registry()) {
        if (t.id.compare(0, filter.size(), filter) != 0) continue;
        ++run;
        std::string msg;
        sc::Timer timer;
        try {
            t.fn();
        } catch (const sc::CheckFailure& e) {
            msg = e.what();
        } catch (const std::exception& e) {
            msg = std::string("exception: ") + e.what();
        }
        double ms = timer.ms();
        bool ok = msg.empty();
        if (!ok) ++failed;
        if (ok)
            std::printf("PASS  %-6s %s (%.1f ms)\n", t.id.c_str(), t.name.c_str(), ms);
        else
            std::printf("FAIL  %-6s %s\n      %s\n", t.id.c_str(), t.name.c_str(), msg.c_str());
        std::fflush(stdout);
        if (csv) csv->row(t.id, t.name, ok ? "PASS" : "FAIL", ms, msg);
    }
    std::error_code ec;
    std::filesystem::remove_all("tests_tmp", ec);
    if (csv) std::printf("Every test with its status and time: %s\n", csv_path.c_str());

    if (run == 0) {
        std::printf("NO TESTS MATCH \"%s\"\n", filter.c_str());
        return 1;
    }
    if (failed) {
        std::printf("%d OF %d TESTS FAILED\n", failed, run);
        return 1;
    }
    std::printf("ALL %d TESTS PASSED\n", run);
    return 0;
}
