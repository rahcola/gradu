class DeBruijn {
public:
  BidirectionalBWTIndex index;
  sdsl::bit_vector first;

  DeBruijn(BidirectionalBWTIndex&&);
};
