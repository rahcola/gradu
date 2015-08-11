#include <iostream>
#include <tuple>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "BidirectionalBWTIndex.hpp"
#include "DeBruijn.hpp"

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

int main(int argc, char *argv[]) {
  if (argc < 3) { return 1; }
  BidirectionalBWTIndex::index_type forward;
  BidirectionalBWTIndex::index_type backward;
  sdsl::construct(forward, argv[1], 1);
  sdsl::construct(backward, argv[2], 1);
  BidirectionalBWTIndex index(std::move(forward), std::move(backward));
  DeBruijn graph(std::move(index));
  std::vector<DeBruijn::size_type> offsets {graph.index.forward.size() / 4,
                                            2 * graph.index.forward.size() / 4,
                                            3 * graph.index.forward.size() / 4,
                                            graph.index.forward.size()};
  DeBruijn::coloring coloring = graph.color(offsets);
  std::cout << "index: "
            << (sdsl::size_in_mega_bytes(graph.index.forward) +
                sdsl::size_in_mega_bytes(graph.index.backward))
            << std::endl;
  std::cout << "first: "
            << (sdsl::size_in_mega_bytes(graph.first) +
                sdsl::size_in_mega_bytes(graph.first_rank) +
                sdsl::size_in_mega_bytes(graph.first_select))
            << std::endl;
  double size = 0;
  for (unsigned int i = 0; i < coloring.size(); i++) {
    size += sdsl::size_in_mega_bytes(std::get<0>(coloring[i]));
    size += sdsl::size_in_mega_bytes(std::get<1>(coloring[i]));
  }
  std::cout << "coloring: " << size << std::endl;

  // for (unsigned int i = 0; i < graph.first.size(); i++) {
  //   std::cout << std::setw(2)
  //             << i
  //             << " "
  //             << graph.first[i]
  //             << " "
  //             << sdsl::extract(graph.index.forward,
  //                              graph.index.forward[i],
  //                              graph.index.forward.size() - 1)
  //             << std::endl;
  // }

  // std::cout << graph.getArc(std::make_tuple(10, 13), 'A')
  //           << " "
  //           << graph.followArc(graph.getArc(std::make_tuple(10, 13), 'A'))
  //           << " "
  //           << graph.getFreq(graph.followArc(graph.getArc(std::make_tuple(10, 13), 'A')))
  //           << std::endl;
  // std::cout << graph.getFreq(graph.followArc(graph.getArc(std::make_tuple(10, 13), 'A')),
  //                            coloring)[0]
  //           << " "
  //           << graph.getFreq(graph.followArc(graph.getArc(std::make_tuple(10, 13), 'A')),
  //                            coloring)[1]
  //           << std::endl;

  return 0;
}
