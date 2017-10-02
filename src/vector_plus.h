#include <functional>
#include <vector>

template <class T>
class VectorPlus
    : public std::
          binary_function<std::vector<T>, std::vector<T>, std::vector<T>> {
 public:
  std::vector<T> operator()(
      const std::vector<T>& lhs, const std::vector<T>& rhs) const {
    std::vector<T> v(lhs.size());
    std::transform(
        lhs.begin(), lhs.end(), rhs.begin(), v.begin(), std::plus<T>());
    return v;
  }
};