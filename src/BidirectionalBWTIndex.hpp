class InternalNodeIterable;

class BidirectionalBWTIndex {
public:
  typedef sdsl::csa_wt<sdsl::wt_hutu<>, 32, 32> index_type;
  typedef index_type::wavelet_tree_type::value_type value_type;
  typedef index_type::wavelet_tree_type::size_type size_type;
  typedef std::tuple<size_type, size_type> interval;
  index_type forward;
  index_type backward;

  BidirectionalBWTIndex(index_type&&, index_type&&);

  std::vector<value_type>&
  enumerateLeft(interval, std::vector<value_type>&) const;
  std::vector<value_type>&
  enumerateRight(interval, std::vector<value_type>&) const;

  std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>
  extendLeft(value_type, interval, interval) const;
  std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>
  extendRight(value_type, interval, interval) const;

  std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>>&
  extendLeftAll(interval, interval, std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>>&) const;
  std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>>&
  extendRightAll(interval, interval, std::vector<std::tuple<BidirectionalBWTIndex::interval, BidirectionalBWTIndex::interval>>&) const;

  bool rightMaximal(interval) const;

  InternalNodeIterable internalNodeIterable();
private:
  mutable std::vector<value_type> symbol_buffer;
  mutable std::vector<size_type> rank_i_buffer;
  mutable std::vector<size_type> rank_j_buffer;
};

class InternalNodeIterator
  : public std::iterator<std::input_iterator_tag,
                         std::tuple<BidirectionalBWTIndex::interval,
                                    unsigned int>> {
public:
  InternalNodeIterator(const BidirectionalBWTIndex &);
  InternalNodeIterator(const BidirectionalBWTIndex &, bool);
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
  std::vector<std::tuple<BidirectionalBWTIndex::interval,
                         BidirectionalBWTIndex::interval>> intervals;
};

class InternalNodeIterable {
public:
  InternalNodeIterable(const BidirectionalBWTIndex &);
  InternalNodeIterator begin();
  InternalNodeIterator end();
private:
  const BidirectionalBWTIndex &index;
};
