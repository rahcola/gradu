class DeBruijn {
public:
  typedef BidirectionalBWTIndex::interval node;
  typedef BidirectionalBWTIndex::interval edge;
  typedef BidirectionalBWTIndex::value_type symbol;
  typedef BidirectionalBWTIndex::size_type size_type;
  typedef sdsl::rrr_vector<> first_vector;
  typedef sdsl::rrr_vector<> coloring_vector;
  typedef std::vector<std::tuple<coloring_vector, coloring_vector::select_1_type>> coloring;

  BidirectionalBWTIndex index;
  first_vector first;

  DeBruijn(BidirectionalBWTIndex&&, unsigned int);

  edge getArc(node, symbol);

  node followArc(edge);

  size_type getFreq(node);
  std::vector<size_type> getFreq(node, coloring&);

  coloring color(std::vector<size_type>&);

  first_vector::rank_1_type first_rank;
  first_vector::select_1_type first_select;
};
