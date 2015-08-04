class DeBruijn {
public:
  typedef BidirectionalBWTIndex::interval node;
  typedef BidirectionalBWTIndex::interval edge;
  typedef BidirectionalBWTIndex::value_type symbol;
  typedef BidirectionalBWTIndex::size_type size_type;

  BidirectionalBWTIndex index;
  sdsl::bit_vector first;

  DeBruijn(BidirectionalBWTIndex&&);

  edge getArc(node, symbol);

  node followArc(node, symbol);

  size_type getFreq(node);
private:
  typedef sdsl::rank_support_v<> rank_support;
  typedef sdsl::select_support_mcl<> select_support;
  rank_support first_rank;
  select_support first_select;
};
