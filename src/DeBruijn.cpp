#include <vector>
#include <sdsl/suffix_arrays.hpp>
#include "BidirectionalBWTIndex.hpp"
#include "DeBruijn.hpp"

DeBruijn::DeBruijn(BidirectionalBWTIndex &&_index)
  : index(std::move(_index)), first(_index.forward.size(), 1) {
  typedef BidirectionalBWTIndex::size_type size_type;
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
}
