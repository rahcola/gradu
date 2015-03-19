#include <iostream>
#include <iterator>
#include <tuple>
#include <stack>
#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

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

class InternalNodeIterator;

class BidirectionalBWTIndex {
public:
  typedef sdsl::csa_wt<sdsl::wt_hutu<sdsl::rrr_vector<127>>, 512, 1024> index_type;
  typedef index_type::wavelet_tree_type::value_type value_type;
  typedef index_type::wavelet_tree_type::size_type size_type;
  typedef std::tuple<size_type, size_type> interval;
  index_type forward;
  index_type backward;

  BidirectionalBWTIndex();
  BidirectionalBWTIndex(std::string, std::string);

  std::vector<value_type>& enumerateLeft(interval, std::vector<value_type>&) const;
  std::vector<value_type>& enumerateRight(interval, std::vector<value_type>&) const;

  std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>
  extendLeft(value_type, interval, interval) const;
  std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>
  extendRight(value_type, interval, interval) const;

  bool rightMaximal(interval) const;

  InternalNodeIterator begin();
  InternalNodeIterator end();
private:
  mutable std::vector<value_type> symbol_buffer;
  mutable std::vector<size_type> rank_buffer;
};

class InternalNodeIterator
  : public std::iterator<std::input_iterator_tag,
                         std::tuple<BidirectionalBWTIndex::interval, unsigned int>> {
public:
  InternalNodeIterator();
  InternalNodeIterator(const BidirectionalBWTIndex &);
  InternalNodeIterator(const InternalNodeIterator &);
  std::tuple<BidirectionalBWTIndex::interval,
             BidirectionalBWTIndex::interval,
             unsigned int> operator*();
  bool operator==(const InternalNodeIterator &);
  bool operator!=(const InternalNodeIterator &);
  InternalNodeIterator &operator++();
  InternalNodeIterator operator++(int);
private:
  const BidirectionalBWTIndex &index;
  std::stack<std::tuple<BidirectionalBWTIndex::interval,
                        BidirectionalBWTIndex::interval,
                        unsigned int>> stack;
  std::vector<BidirectionalBWTIndex::value_type> symbols;
};

BidirectionalBWTIndex::BidirectionalBWTIndex() { }

BidirectionalBWTIndex::BidirectionalBWTIndex(std::string fforward,
                                             std::string fbackward) {
  sdsl::construct(forward, fforward, 1);
  sdsl::construct(backward, fbackward, 1);
  symbol_buffer.resize(forward.sigma);
  rank_buffer.resize(forward.sigma);
}

std::vector<BidirectionalBWTIndex::value_type>&
BidirectionalBWTIndex::enumerateLeft(interval ij,
                                     std::vector<value_type> &out) const {
  out.resize(forward.sigma);
  size_type i, j, k;
  std::tie(i, j) = ij;
  sdsl::interval_symbols(forward.wavelet_tree, i, j + 1,
                         k, out, rank_buffer, rank_buffer);
  out.resize(k);
  return out;
}

std::vector<BidirectionalBWTIndex::value_type>&
BidirectionalBWTIndex::enumerateRight(interval ij,
                                      std::vector<value_type> &out) const {
  out.resize(backward.sigma);
  size_type i, j, k;
  std::tie(i, j) = ij;
  sdsl::interval_symbols(backward.wavelet_tree, i, j + 1,
                         k, out, rank_buffer, rank_buffer);
  out.resize(k);
  return out;
}

std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>
BidirectionalBWTIndex::extendLeft(value_type c, interval ij, interval pq) const {
  size_type i, j, ii, jj;
  size_type p = std::get<0>(pq);
  std::tie(i, j) = ij;
  sdsl::backward_search(forward, i, j, c, ii, jj);
  size_type k = std::get<1>(forward.wavelet_tree.lex_count(i, j + 1, c));
  return std::make_tuple(std::make_tuple(ii, jj),
                         std::make_tuple(p + k, p + k + jj - ii));
}

std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>
BidirectionalBWTIndex::extendRight(value_type c, interval ij, interval pq) const {
  size_type p, q, pp, qq;
  size_type i = std::get<0>(ij);
  std::tie(p, q) = pq;
  sdsl::backward_search(backward, p, q, c, pp, qq);
  size_type k = std::get<1>(backward.wavelet_tree.lex_count(p, q + 1, c));
  return std::make_tuple(std::make_tuple(i + k, i + k + qq - pp),
                         std::make_tuple(pp, qq));
}

bool BidirectionalBWTIndex::rightMaximal(interval ij) const {
  size_type i, j, k;
  std::tie(i, j) = ij;
  sdsl::interval_symbols(backward.wavelet_tree, i, j + 1,
                         k, symbol_buffer, rank_buffer, rank_buffer);
  return k > 1;
}

InternalNodeIterator BidirectionalBWTIndex::begin() {
  return InternalNodeIterator(*this);
}

InternalNodeIterator BidirectionalBWTIndex::end() {
  return InternalNodeIterator();
}

InternalNodeIterator::InternalNodeIterator()
  : index(BidirectionalBWTIndex()) { }

InternalNodeIterator::InternalNodeIterator(const BidirectionalBWTIndex &_index)
  : index(_index) {
  stack.push(std::make_tuple(std::make_tuple(0, index.forward.size() - 1),
                             std::make_tuple(0, index.backward.size() - 1),
                             0));
}

InternalNodeIterator::InternalNodeIterator(const InternalNodeIterator &other)
  : index(other.index), stack(other.stack) {
}

bool InternalNodeIterator::operator==(const InternalNodeIterator &other) {
  return stack == other.stack;
}

bool InternalNodeIterator::operator!=(const InternalNodeIterator &other) {
  return !(*this == other);
}

std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval,
           unsigned int>
InternalNodeIterator::operator*() {
  return stack.top();
}

InternalNodeIterator &InternalNodeIterator::operator++() {
  typedef BidirectionalBWTIndex::interval interval;
  interval ij, pq;
  BidirectionalBWTIndex::size_type d;
  std::tie(ij, pq, d) = stack.top();
  stack.pop();

  for (auto c: index.enumerateLeft(ij, symbols)) {
    interval iijj, ppqq;
    std::tie(iijj, ppqq) = index.extendLeft(c, ij, pq);
    if (index.rightMaximal(ppqq)) {
      stack.push(std::make_tuple(iijj, ppqq, d + 1));
    }
  }

  return *this;
}

InternalNodeIterator InternalNodeIterator::operator++(int dummy) {
  InternalNodeIterator it(*this);
  this->operator++();
  return it;
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

int main(int argc, char *argv[]) {
  typedef BidirectionalBWTIndex::interval interval;
  typedef BidirectionalBWTIndex::value_type value_type;
  typedef BidirectionalBWTIndex::size_type size_type;

  if (argc < 3) { return 1; }
  std::string fforward(argv[1]);
  std::string fbackward(argv[2]);
  BidirectionalBWTIndex index(fforward, fbackward);
  std::vector<value_type> symbols;
  sdsl::bit_vector first(index.forward.size(), 1);

  for(auto it : index) {
    interval ij, pq;
    unsigned int d;
    std::tie(ij, pq, d) = it;
    if (d >= 2) {
      bool f = true;
      for (auto c: index.enumerateRight(pq, symbols)) {
        if (f) { f = false; continue; }
        interval iijj = std::get<0>(index.extendRight(c, ij, pq));
        size_type ii = std::get<0>(iijj);
        first[ii] = first[ii] ^ 1;
      }
    }
  }

  sdsl::rank_support_v5<> first_rank(&first);
  sdsl::select_support_mcl<> first_select(&first);
  std::cout << followArc(index, first, first_rank, first_select,
                         std::make_tuple(6, 8), 'G')
            << std::endl;
  return 0;
}
