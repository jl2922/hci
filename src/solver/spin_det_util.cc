#include "spin_det_util.h"

// Converting from HF diffs representation to occupied orbitals.
std::vector<int> SpinDetUtil::get_occupied_orbitals(
    const data::SpinDeterminant& spin_det) {
  std::vector<int> occ;
  const int n_hf_elecs = spin_det.n_hf_elecs();
  occ.reserve(n_hf_elecs);
  int level = 0;
  for (const int v_orb : spin_det.v_holes()) {
    while (level < v_orb) {
      occ.push_back(level);
      level++;
    }
    level = v_orb + 1;
  }
  while (level < n_hf_elecs) {
    occ.push_back(level);
    level++;
  }
  for (const int c_orb : spin_det.c_elecs()) {
    occ.push_back(c_orb);
  }
  return occ;
};

bool SpinDetUtil::is_occupied(
    const data::SpinDeterminant& spin_det, const int orb) {
  const int n_hf_elecs = spin_det.n_hf_elecs();
  if (orb < n_hf_elecs) {
    for (const int v_orb : spin_det.v_holes()) {
      if (v_orb == orb) return false;
      if (v_orb > orb) return true;
    }
    return true;
  } else {
    for (const int c_orb : spin_det.c_elecs()) {
      if (c_orb == orb) return true;
      if (c_orb > orb) return false;
    }
    return false;
  }
}

int SpinDetUtil::get_n_lower_elecs(
    const data::SpinDeterminant& spin_det, const int orb) {
  const int n_hf_elecs = spin_det.n_hf_elecs();
  if (orb < n_hf_elecs) {
    int res = orb;
    for (const int v_orb : spin_det.v_holes()) {
      if (v_orb == orb) return 0;
      if (v_orb > orb) return res;
      res--;
    }
    return res;
  } else {
    int res = n_hf_elecs - spin_det.v_holes_size();
    for (const int c_orb : spin_det.c_elecs()) {
      if (c_orb == orb) return res;
      if (c_orb > orb) return 0;
      res++;
    }
    return 0;
  }
}

int SpinDetUtil::get_highest_orbital(const data::SpinDeterminant& spin_det) {
  const int n_hf_elecs = spin_det.n_hf_elecs();
  const int n_c_elecs = spin_det.c_elecs_size();
  const int n_v_holes = spin_det.v_holes_size();
  if (n_c_elecs > 0) {
    return spin_det.c_elecs(n_c_elecs - 1);
  } else if (n_v_holes == 0) {
    return n_hf_elecs - 1;
  }
  throw std::runtime_error("n elecs not conserved.");
}

void SpinDetUtil::set_occupation(
    data::SpinDeterminant* spin_det, const int orb, const bool occ) {
  const int n_hf_elecs = spin_det->n_hf_elecs();
  if (orb < n_hf_elecs) {
    if (occ) {
      remove_orbital(spin_det->mutable_v_holes(), orb);
    } else {
      insert_orbital(spin_det->mutable_v_holes(), orb);
    }
  } else {
    if (occ) {
      insert_orbital(spin_det->mutable_c_elecs(), orb);
    } else {
      remove_orbital(spin_det->mutable_c_elecs(), orb);
    }
  }
}

std::vector<int> SpinDetUtil::get_eor(
    const data::SpinDeterminant& lhs, const data::SpinDeterminant& rhs) {
  std::vector<int> eor_orbs;
  eor_orbs.reserve(5);

  // Process v holes.
  int lhs_ptr = 0;
  int rhs_ptr = 0;
  int lhs_size = lhs.v_holes_size();
  int rhs_size = rhs.v_holes_size();
  while (lhs_ptr < lhs_size && rhs_ptr < rhs_size) {
    if (lhs.v_holes(lhs_ptr) == rhs.v_holes(rhs_ptr)) {
      lhs_ptr++;
      rhs_ptr++;
    } else if (lhs.v_holes(lhs_ptr) > rhs.v_holes(rhs_ptr)) {
      eor_orbs.push_back(rhs.v_holes(rhs_ptr));
      rhs_ptr++;
    } else {
      eor_orbs.push_back(lhs.v_holes(lhs_ptr));
      lhs_ptr++;
    }
    if (eor_orbs.size() > 4) return eor_orbs;
  }
  while (lhs_ptr < lhs_size) {
    eor_orbs.push_back(lhs.v_holes(lhs_ptr));
    if (eor_orbs.size() > 4) return eor_orbs;
    lhs_ptr++;
  }
  while (rhs_ptr < rhs_size) {
    eor_orbs.push_back(rhs.v_holes(rhs_ptr));
    if (eor_orbs.size() > 4) return eor_orbs;
    rhs_ptr++;
  }

  // Process c elecs.
  lhs_ptr = 0;
  rhs_ptr = 0;
  lhs_size = lhs.c_elecs_size();
  rhs_size = rhs.c_elecs_size();
  while (lhs_ptr < lhs_size && rhs_ptr < rhs_size) {
    if (lhs.c_elecs(lhs_ptr) == rhs.c_elecs(rhs_ptr)) {
      lhs_ptr++;
      rhs_ptr++;
    } else if (lhs.c_elecs(lhs_ptr) > rhs.c_elecs(rhs_ptr)) {
      eor_orbs.push_back(rhs.c_elecs(rhs_ptr));
      rhs_ptr++;
    } else {
      eor_orbs.push_back(lhs.c_elecs(lhs_ptr));
      lhs_ptr++;
    }
    if (eor_orbs.size() > 4) return eor_orbs;
  }
  while (lhs_ptr < lhs_size) {
    eor_orbs.push_back(lhs.c_elecs(lhs_ptr));
    if (eor_orbs.size() > 4) return eor_orbs;
    lhs_ptr++;
  }
  while (rhs_ptr < rhs_size) {
    eor_orbs.push_back(rhs.c_elecs(rhs_ptr));
    if (eor_orbs.size() > 4) return eor_orbs;
    rhs_ptr++;
  }

  return eor_orbs;
}

void SpinDetUtil::insert_orbital(
    google::protobuf::RepeatedField<google::protobuf::int32>* t_orbs,
    const int orb) {
  t_orbs->Add();
  const int n_t_orbs = t_orbs->size();
  int i;
  for (i = n_t_orbs - 1; i >= 1; i--) {
    const int t_orb = t_orbs->Get(i - 1);
    if (t_orb == orb) throw std::runtime_error("orbital already occupied");
    if (t_orb < orb) break;
    *(t_orbs->Mutable(i)) = t_orbs->Get(i - 1);
  }
  *(t_orbs->Mutable(i)) = orb;
}

void SpinDetUtil::remove_orbital(
    google::protobuf::RepeatedField<google::protobuf::int32>* t_orbs,
    const int orb) {
  bool found = false;
  const int n_t_orbs = t_orbs->size();
  for (int i = 0; i < n_t_orbs; i++) {
    const int t_orb = t_orbs->Get(i);
    if (found) *(t_orbs->Mutable(i - 1)) = t_orb;
    if (t_orb == orb) found = true;
    if (t_orb > orb && !found) break;
  }
  if (found) t_orbs->RemoveLast();
}
