#pragma once

#include <stdexcept>


namespace fgt {


class fgt_error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};


class dimension_mismatch : public fgt_error {
public:
    using fgt_error::fgt_error;
};
}
