#ifndef SPIN_DET_UTIL_H_
#define SPIN_DET_UTIL_H_

#include <vector>
#include "../data.pb.h"

class SpinDetUtil {
 public:
  static std::vector<int> get_occupied_orbitals(
      const data::SpinDeterminant& spin_det);

  static bool is_occupied(const data::SpinDeterminant& spin_det, const int orb);

  static int get_n_lower_elecs(
      const data::SpinDeterminant& spin_det, const int orb);

  static int get_highest_orbital(const data::SpinDeterminant& spin_det);

  static void set_occupation(
      data::SpinDeterminant* spin_det, const int orb, const bool occ);

  static std::vector<int> get_eor(
      const data::SpinDeterminant& lhs, const data::SpinDeterminant& rhs);

 private:
  static void insert_orbital(
      google::protobuf::RepeatedField<google::protobuf::int32>* orbs,
      const int orb);

  static void remove_orbital(
      google::protobuf::RepeatedField<google::protobuf::int32>* orbs,
      const int orb);
};

#endif
