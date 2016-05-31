#ifndef CAFFE_UNIFIED_DATA_READER_HPP_
#define CAFFE_UNIFIED_DATA_READER_HPP_

#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include "caffe/common.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/util/blocking_queue.hpp"
#include "caffe/image_db.hpp"
#include "caffe/proto/caffe.pb.h"

namespace caffe {

// Unified data input is a tsv file, each column map to an instance
struct UnifiedSample {
  vector<ImageData>* image_data;
  vector<string> data;
  float label;
};

/**
 * @brief Reads data from a source to queues available to data layers.
 * A single reading thread is created per source, even if multiple solvers
 * are running in parallel, e.g. for multi-GPU training. This makes sure
 * databases are read sequentially, and that each solver accesses a different
 * subset of the database. Data is distributed to solvers in a round-robin
 * way to keep parallel training deterministic.
 * The data source is a tsv file, each column maps to an output of this layer.
 */
class UnifiedDataReader {
 public:
  explicit UnifiedDataReader(const LayerParameter& param);
  ~UnifiedDataReader();

  inline BlockingQueue<UnifiedSample*>& free() const {
    return queue_pair_->free_;
  }
  inline BlockingQueue<UnifiedSample*>& full() const {
    return queue_pair_->full_;
  }

 protected:
  // Queue pairs are shared between a body and its readers
  class QueuePair {
   public:
    explicit QueuePair(int size);
    ~QueuePair();

    BlockingQueue<UnifiedSample*> free_;
    BlockingQueue<UnifiedSample*> full_;

  DISABLE_COPY_AND_ASSIGN(QueuePair);
  };

  // A single body is created per source
  class Body : public InternalThread {
   public:
    explicit Body(const LayerParameter& param);
    virtual ~Body();

   protected:
    void InternalThreadEntry();
    void read_one(QueuePair* qp); 
    void PreloadIUB(UnifiedSample* pstr);
    void PrefetchOneInMemory(UnifiedSample* pstr);
    void PrefetchOneFromReadingFile(UnifiedSample* pstr);
    void ShuffleSamples();

    const LayerParameter param_;
    BlockingQueue<shared_ptr<QueuePair> > new_queue_pairs_;

    std::ifstream data_stream_;
    string source_;
    bool load_in_memory_;
    bool shuffle_;
    int sample_id_;
    int batch_size_;
    shared_ptr<Caffe::RNG> prefetch_rng_;

    // Store the whole set for load_in_memory setting
    vector<shared_ptr<UnifiedSample>> samples_;
    vector<shared_ptr<IImageDB>> imageDBs;

    friend class UnifiedDataReader;

  DISABLE_COPY_AND_ASSIGN(Body);
  };

  // Category_Based_Sample_Body is to get k samples from each class
  class Category_Based_Sample_Body : public Body {
  public:
    explicit Category_Based_Sample_Body(const LayerParameter& param);

  protected:
    void InternalThreadEntry();
    //void read_one(QueuePair* qp);

    int samples_per_class_;
    shared_ptr<Caffe::RNG> sample_rng_;

    // maintain the samples in each class
    std::unordered_map<float, vector<int>> class_sample_list_;
    // number of samples in each class, this is for weighted sampling
    vector<int> class_sample_list_len_;
    // log the erased classes to evaluate class
    vector<float> erased_class_;
    // sample_idx_in_ori_list_[i] is the indexes of that sample in the list before erasing
    vector<int> sample_idx_in_ori_list_;
    // give each label an index to handle float or non-continous label
    std::unordered_map<int, float> idx_to_class_;

    friend class UnifiedDataReader;

    DISABLE_COPY_AND_ASSIGN(Category_Based_Sample_Body);
  };

  // A source is uniquely identified by its layer name + path, in case
  // the same database is read from two different locations in the net.
  static inline string source_key(const LayerParameter& param) {
    return param.name() + ":" + param.unified_data_param().source();
  }

  const shared_ptr<QueuePair> queue_pair_;
  shared_ptr<Body> body_;
  static map<const string, boost::weak_ptr<UnifiedDataReader::Body> > bodies_;

DISABLE_COPY_AND_ASSIGN(UnifiedDataReader);
};

}  // namespace caffe

#endif  // CAFFE_UNIFIED_DATA_READER_HPP_
