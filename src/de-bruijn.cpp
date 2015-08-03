#include <iostream>
#include <tuple>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "BidirectionalBWTIndex.hpp"

namespace aux{
template<std::size_t...> struct seq{};

template<std::size_t N, std::size_t... Is>
struct gen_seq : gen_seq<N-1, N-1, Is...>{};

template<std::size_t... Is>
struct gen_seq<0, Is...> : seq<Is...>{};

template<class Ch, class Tr, class Tuple, std::size_t... Is>
void print_tuple(std::basic_ostream<Ch,Tr>& os, Tuple const& t, seq<Is...>){
  using swallow = int[];
  (void)swallow{0, (void(os << (Is == 0? "" : ", ") << std::get<Is>(t)), 0)...};
}
} // aux::

template<class Ch, class Tr, class... Args>
auto operator<<(std::basic_ostream<Ch, Tr>& os, std::tuple<Args...> const& t)
    -> std::basic_ostream<Ch, Tr>&
{
  os << "(";
  aux::print_tuple(os, t, aux::gen_seq<sizeof...(Args)>());
  return os << ")";
}

BidirectionalBWTIndex::interval
getArc(BidirectionalBWTIndex &index,
       BidirectionalBWTIndex::interval v,
       BidirectionalBWTIndex::value_type c) {
  BidirectionalBWTIndex::size_type i, j, ii, jj;
  std::tie(i, j) = v;
  sdsl::backward_search(index.forward, i, j, c, ii, jj);
  return std::make_tuple(ii, jj);
}

BidirectionalBWTIndex::interval
followArc(BidirectionalBWTIndex &index,
          sdsl::bit_vector &first,
          sdsl::rank_support &first_rank,
          sdsl::select_support &first_select,
          BidirectionalBWTIndex::interval v,
          BidirectionalBWTIndex::value_type c) {
  BidirectionalBWTIndex::size_type i, j;
  std::tie(i, j) = getArc(index, v, c);
  BidirectionalBWTIndex::size_type ii = i;
  BidirectionalBWTIndex::size_type jj = j;
  if (!first[i]) {
    ii = first_select.select(first_rank.rank(i));
  }
  if (j < first.size() - 1) {
    jj = first_select.select(first_rank.rank(j) + 1);
  }
  return std::make_tuple(ii, jj);
}

std::vector<int>
getFreqs(BidirectionalBWTIndex &index,
         std::vector<BidirectionalBWTIndex::value_type> color_borders,
         BidirectionalBWTIndex::interval v) {
  std::vector<int> freqs(color_borders.size());
  BidirectionalBWTIndex::size_type i, j;
  std::tie(i, j) = v;
  for (; i < j; i++) {
    for (auto border = color_borders.begin();
         border != color_borders.end();
         border++) {
      if (index.forward[i] < *border) {
        freqs[border - color_borders.begin()] += 1;
        break;
      }
    }
  }
  return freqs;
}

int main(int argc, char *argv[]) {
  typedef BidirectionalBWTIndex::interval interval;
  typedef BidirectionalBWTIndex::size_type size_type;

  if (argc < 3) { return 1; }
  std::string fforward(argv[1]);
  std::string fbackward(argv[2]);
  // BidirectionalBWTIndex::index_type f;
  // BidirectionalBWTIndex::index_type b;
  // sdsl::load_from_file(f, fforward);
  // sdsl::load_from_file(b, fbackward);
  // BidirectionalBWTIndex index(f, b);
  BidirectionalBWTIndex index(fforward, fbackward);
  std::vector<std::tuple<interval, interval>> intervals;
  sdsl::bit_vector first(index.forward.size(), 1);

  for(auto it : index.internalNodeIterable()) {
    interval ij, pq;
    unsigned int d;
    std::tie(ij, pq, d) = it;
    if (d >= 2) {
      bool first_c = true;
      for (auto is : index.extendRightAll(ij, pq, intervals)) {
        if (first_c) { first_c = false; continue; }
        size_type i = std::get<0>(std::get<0>(is));
        first[i] = first[i] ^ 1;
      }
    }
  }

  for (unsigned int i = 0; i < first.size(); i++) {
    std::cout << first[i]
              << " "
              << sdsl::extract(index.forward,
                               index.forward[i], index.forward.size() - 1)
              << std::endl;
  }

  return 0;
}
