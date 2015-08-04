#include <vector>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/util.hpp>
#include "BidirectionalBWTIndex.hpp"
#include "DeBruijn.hpp"

DeBruijn::DeBruijn(BidirectionalBWTIndex &&_index)
  : index(std::move(_index)), first(_index.forward.size(), 1) {
  typedef DeBruijn::size_type size_type;
  typedef BidirectionalBWTIndex::interval interval;
  std::vector<std::tuple<interval, interval>> intervals;
  for(auto it : index.internalNodeIterable()) {
    interval ij, pq;
    unsigned int d;
    std::tie(ij, pq, d) = it;
    if (d >= 2) {
      auto it = index.extendRightAll(ij, pq, intervals);
      std::for_each(it.begin() + 1, it.end(),
                    [this](std::tuple<interval, interval> &t) {
                      size_type i = std::get<0>(std::get<0>(t));
                      first[i] = first[i] ^ 1;
                    });
    }
  }
  sdsl::util::init_support(first_rank, &first);
  sdsl::util::init_support(first_select, &first);
}

DeBruijn::edge DeBruijn::getArc(DeBruijn::node v, DeBruijn::symbol c) {
  DeBruijn::size_type i, j, ii, jj;
  std::tie(i, j) = v;
  sdsl::backward_search(index.forward, i, j, c, ii, jj);
  return std::make_tuple(ii, jj);
}

DeBruijn::node DeBruijn::followArc(DeBruijn::node v, DeBruijn::symbol c) {
  DeBruijn::size_type i, j;
  std::tie(i, j) = getArc(v, c);
  DeBruijn::size_type ii = i;
  DeBruijn::size_type jj = j;
  if (!first[i]) {
    ii = first_select.select(first_rank.rank(i));
  }
  if (j < first.size() - 1) {
    jj = first_select.select(first_rank.rank(j) + 1);
  }
  return std::make_tuple(ii, jj);
}

DeBruijn::size_type DeBruijn::getFreq(DeBruijn::node v) {
  return std::get<1>(v) - std::get<0>(v);
}
