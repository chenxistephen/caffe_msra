#ifndef CAFFE_RNG_CPP_HPP_
#define CAFFE_RNG_CPP_HPP_

#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <unordered_set>

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int.hpp"

#include "caffe/common.hpp"

namespace caffe {

typedef boost::mt19937 rng_t;

inline rng_t* caffe_rng() {
  return static_cast<caffe::rng_t*>(Caffe::rng_stream().generator());
}

// Fisher–Yates algorithm
template <class RandomAccessIterator, class RandomGenerator>
inline void shuffle(RandomAccessIterator begin, RandomAccessIterator end,
                    RandomGenerator* gen) {
  typedef typename std::iterator_traits<RandomAccessIterator>::difference_type
      difference_type;
  typedef typename boost::uniform_int<difference_type> dist_type;

  difference_type length = std::distance(begin, end);
  if (length <= 0) return;

  for (difference_type i = length - 1; i > 0; --i) {
    dist_type dist(0, i);
    std::iter_swap(begin + i, begin + dist(*gen));
  }
}

template <class RandomAccessIterator>
inline void shuffle(RandomAccessIterator begin, RandomAccessIterator end) {
  shuffle(begin, end, caffe_rng());
}

// return a random permutation of a sequence. If not weighted, each element is equally occured at different places; If weighted, the more weighted index have the chance to occur earlier;
class RandomPermutation {
private:
  int sample_range_;
  std::unordered_map<int, int> reorder_;
  vector<int> sample_index_;
  std::unordered_set<int> used_;
  inline int ReorderedWithMap(int id)
  {
    auto it = reorder_.find(id);
    return (it != reorder_.end()) ? it->second : id;
  }
public:
  RandomPermutation(int range) : sample_range_(range) {
    for (int i = 0; i <= range; i++)
      sample_index_.push_back(i);
    used_.clear();
  }
  RandomPermutation(const vector<int>& index_weight) {
    for (int i = 0; i < index_weight.size(); i++) {
      for (int j = 0; j < index_weight[i]; j++) {
        sample_index_.push_back(i);
      }
    }
    sample_range_ = sample_index_.size() - 1;
    used_.clear();
  }
  inline void Reset() { reorder_.clear(); used_.clear(); }
  int GetNext(caffe::rng_t* rng) {
    int rid = -1, id = -1;
    while (true) {
      rid = boost::uniform_int<int>(0, sample_range_)(*rng);
      id = ReorderedWithMap(rid);
      reorder_[rid] = ReorderedWithMap(sample_range_--);
      if (used_.insert(sample_index_[id]).second)
        break;
    }
    CHECK(id != -1) << "should get a valid id.";
    return sample_index_[id];
  }
};


}  // namespace caffe

#endif  // CAFFE_RNG_HPP_
