#pragma once

#include <iostream>
#include <string>

namespace meow::test {

inline int failures = 0;

inline void check(bool condition, const std::string& message) {
    if (condition) {
        std::cout << "PASS: " << message << '\n';
    } else {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

inline int summary() {
    if (failures != 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    return 0;
}

}  // namespace meow::test
