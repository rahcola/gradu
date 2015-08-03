#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "BidirectionalBWTIndex.hpp"

BidirectionalBWTIndex::size_type
count_less(std::vector<BidirectionalBWTIndex::size_type> &,
           std::vector<BidirectionalBWTIndex::size_type> &,
           BidirectionalBWTIndex::size_type);

BidirectionalBWTIndex::BidirectionalBWTIndex() { }

BidirectionalBWTIndex::BidirectionalBWTIndex(std::string fforward,
                                             std::string fbackward) {
  sdsl::construct(forward, fforward, 1);
  sdsl::construct(backward, fbackward, 1);
  symbol_buffer.resize(forward.sigma);
  rank_i_buffer.resize(forward.sigma);
  rank_j_buffer.resize(forward.sigma);
}

BidirectionalBWTIndex::BidirectionalBWTIndex(BidirectionalBWTIndex::index_type _forward,
                                             BidirectionalBWTIndex::index_type _backward) {
  forward = _forward;
  backward = _backward;
  symbol_buffer.resize(forward.sigma);
  rank_i_buffer.resize(forward.sigma);
  rank_j_buffer.resize(forward.sigma);
}

std::vector<BidirectionalBWTIndex::value_type>&
BidirectionalBWTIndex::enumerateLeft(interval ij,
                                     std::vector<value_type> &out) const {
  out.resize(forward.sigma);
  size_type i, j, k;
  std::tie(i, j) = ij;
  sdsl::interval_symbols(forward.wavelet_tree, i, j + 1,
                         k, out, rank_i_buffer, rank_j_buffer);
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
                         k, out, rank_i_buffer, rank_j_buffer);
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

std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>>&
BidirectionalBWTIndex::extendLeftAll(interval ij, interval pq, std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>> &out) const {
  size_type i, j;
  std::tie(i, j) = ij;
  size_type symbol_count = 0;
  size_type p = std::get<0>(pq);
  sdsl::interval_symbols(forward.wavelet_tree, i, j + 1,
                         symbol_count, symbol_buffer,
                         rank_i_buffer, rank_j_buffer);
  out.resize(symbol_count);
  for (std::vector<size_type>::size_type ix = 0; ix < symbol_count; ix++) {
    size_type c = backward.C[backward.char2comp[symbol_buffer[ix]]];
    size_type ii = c + rank_i_buffer[ix];
    size_type jj = c + rank_j_buffer[ix] - 1;
    size_type k = count_less(rank_i_buffer, rank_j_buffer, ix);
    out[ix] = std::make_tuple(std::make_tuple(ii, jj),
                              std::make_tuple(p + k, p + k + jj - ii));
  }
  return out;
}

std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>>&
BidirectionalBWTIndex::extendRightAll(interval ij, interval pq, std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>> &out) const {
  size_type p, q;
  std::tie(p, q) = pq;
  size_type symbol_count = 0;
  size_type i = std::get<0>(ij);
  sdsl::interval_symbols(backward.wavelet_tree, p, q + 1,
                         symbol_count, symbol_buffer,
                         rank_i_buffer, rank_j_buffer);
  out.resize(symbol_count);
  for (std::vector<size_type>::size_type ix = 0; ix < symbol_count; ix++) {
    size_type c = backward.C[backward.char2comp[symbol_buffer[ix]]];
    size_type pp = c + rank_i_buffer[ix];
    size_type qq = c + rank_j_buffer[ix] - 1;
    size_type k = count_less(rank_i_buffer, rank_j_buffer, ix);
    out[ix] = std::make_tuple(std::make_tuple(i + k, i + k + qq - pp),
                              std::make_tuple(pp, qq));
  }
  return out;
}

bool BidirectionalBWTIndex::rightMaximal(interval ij) const {
  size_type i, j, k;
  std::tie(i, j) = ij;
  sdsl::interval_symbols(backward.wavelet_tree, i, j + 1,
                         k, symbol_buffer, rank_i_buffer, rank_j_buffer);
  return k > 1;
}

BidirectionalBWTIndex::size_type
count_less(std::vector<BidirectionalBWTIndex::size_type> &rank_i_buffer,
           std::vector<BidirectionalBWTIndex::size_type> &rank_j_buffer,
           BidirectionalBWTIndex::size_type ix) {
  BidirectionalBWTIndex::size_type k = 0;
  for (std::vector<BidirectionalBWTIndex::size_type>::size_type jx = 0; jx < ix; jx++) {
    k += rank_j_buffer[jx] - rank_i_buffer[jx];
  }
  return k;
}

InternalNodeIterable BidirectionalBWTIndex::internalNodeIterable() {
  return InternalNodeIterable(*this);
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

  for (auto i : index.extendLeftAll(ij, pq, intervals)) {
    if (index.rightMaximal(std::get<1>(i))) {
      stack.push(std::make_tuple(std::get<0>(i), std::get<1>(i), d + 1));
    }
  }

  return *this;
}

InternalNodeIterator InternalNodeIterator::operator++(int dummy) {
  InternalNodeIterator it(*this);
  this->operator++();
  return it;
}

InternalNodeIterable::InternalNodeIterable(const BidirectionalBWTIndex &_index)
  : index(_index) { }

InternalNodeIterator InternalNodeIterable::begin() {
  return InternalNodeIterator(index);
}

InternalNodeIterator InternalNodeIterable::end() {
  return InternalNodeIterator();
}
