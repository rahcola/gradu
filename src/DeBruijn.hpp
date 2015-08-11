class DeBruijn {
public:
  typedef BidirectionalBWTIndex::interval node;
  typedef BidirectionalBWTIndex::interval edge;
  typedef BidirectionalBWTIndex::value_type symbol;
  typedef BidirectionalBWTIndex::size_type size_type;
  typedef std::vector<std::tuple<sdsl::sd_vector<>, sdsl::select_support_sd<>>> coloring;

  BidirectionalBWTIndex index;
  sdsl::sd_vector<> first;

  DeBruijn(BidirectionalBWTIndex&&);

  edge getArc(node, symbol);

  node followArc(edge);

  size_type getFreq(node);
  std::vector<size_type> getFreq(node, coloring&);

  coloring color(std::vector<size_type>&);

  typedef sdsl::rank_support_sd<> rank_support;
  typedef sdsl::select_support_sd<> select_support;
  rank_support first_rank;
  select_support first_select;
};
