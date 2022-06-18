#pragma once

#include "emp/base/vector.hpp"

namespace selection {

class BaseSelect {
protected:
  emp::vector<size_t> selected;
  std::string name="BaseSelect";

public:
  virtual ~BaseSelect() = default;

  /// Select n individuals (update selected)
  virtual emp::vector<size_t>& operator()(size_t n) = 0;

  emp::vector<size_t>& GetSelected() { return selected; }
  const emp::vector<size_t>& GetSelected() const { return selected; }
};

}