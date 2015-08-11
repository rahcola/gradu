#include <vector>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/util.hpp>
#include "BidirectionalBWTIndex.hpp"
#include "DeBruijn.hpp"

DeBruijn::DeBruijn(BidirectionalBWTIndex &&_index)
  : index(std::move(_index)) {
  typedef DeBruijn::size_type size_type;
  typedef BidirectionalBWTIndex::interval interval;
  sdsl::bit_vector uncompressed_first(index.forward.size(), 1);
  std::vector<std::tuple<interval, interval>> intervals;
  for(auto it : index.internalNodeIterable()) {
    interval ij, pq;
    unsigned int d;
    std::tie(ij, pq, d) = it;
    if (d >= 2) {
      auto it = index.extendRightAll(ij, pq, intervals);
      std::for_each(it.begin() + 1, it.end(),
                    [&uncompressed_first](std::tuple<interval, interval> &t) {
                      size_type i = std::get<0>(std::get<0>(t));
                      uncompressed_first[i] = uncompressed_first[i] ^ 1;
                    });
    }
  }
  sdsl::sd_vector<> tmp(uncompressed_first);
  first.swap(tmp);
  sdsl::util::init_support(first_rank, &first);
  sdsl::util::init_support(first_select, &first);
}

DeBruijn::edge DeBruijn::getArc(DeBruijn::node v, DeBruijn::symbol c) {
  DeBruijn::size_type i, j, ii, jj;
  std::tie(i, j) = v;
  sdsl::backward_search(index.forward, i, j, c, ii, jj);
  return std::make_tuple(ii, jj);
}

DeBruijn::node DeBruijn::followArc(DeBruijn::node v) {
  DeBruijn::size_type i, j;
  std::tie(i, j) = v;
  if (!first[i]) {
    i = first_select.select(first_rank.rank(i));
  }
  if (j < first.size() - 1) {
    j = first_select.select(first_rank.rank(j) + 1);
  }
  return std::make_tuple(i, j);
}

DeBruijn::size_type DeBruijn::getFreq(DeBruijn::node v) {
  return std::get<1>(v) - std::get<0>(v);
}

std::vector<DeBruijn::size_type> DeBruijn::getFreq(DeBruijn::node v,
                                                   DeBruijn::coloring &coloring) {
  size_type start = first_rank.rank(std::get<0>(v)) + 1;
  size_type end = first_rank.rank(std::get<1>(v)) + 1;
  std::vector<size_type> freqs;
  for (auto c : coloring) {
    freqs.push_back((std::get<1>(c).select(end) + 1 - end) -
                    (std::get<1>(c).select(start) + 1 - start));
  }
  return freqs;
}

DeBruijn::coloring DeBruijn::color(std::vector<DeBruijn::size_type> &offsets) {
  size_type m = index.forward.size();
  rank_support::size_type n = first_rank(m) + 1;
  std::vector<sdsl::bit_vector> bits(offsets.size(), sdsl::bit_vector(m + n, 0));
  std::vector<sdsl::bit_vector::size_type> sums(offsets.size(), 0);
  for (size_type i = 0; i < m; i++) {
    if (first[i]) {
      for (unsigned int j = 0; j < bits.size(); j++) {
        bits[j][sums[j]] = 1;
        sums[j]++;
      }
    }
    size_type position = index.forward[i];
    for (unsigned int j = 0; j < offsets.size(); j++) {
      if (position < offsets[j]) {
        sums[j]++;
        break;
      }
    }
  }
  coloring c;
  for (unsigned int i = 0; i < offsets.size(); i++) {
    bits[i][sums[i]] = 1;
    c.emplace_back(sdsl::sd_vector<> (bits[i]), sdsl::select_support_sd<> ());
  }
  for (unsigned int i = 0; i < offsets.size(); i++) {
    sdsl::util::init_support(std::get<1>(c[i]), &std::get<0>(c[i]));
  }
  return c;
}
