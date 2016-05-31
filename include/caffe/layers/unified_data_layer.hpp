#ifndef CAFFE_UNIFIED_DATA_LAYER_HPP_
#define CAFFE_UNIFIED_DATA_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/data_transformer.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/layer.hpp"
#include "caffe/data_loader.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/blocking_queue.hpp"
#include "caffe/unified_data_reader.hpp"
#include "caffe/internal_thread.hpp"

namespace caffe {

template <typename Dtype>
class UnifiedBatch {
public:
  vector<shared_ptr<Blob<Dtype>>> data_;
  Blob<Dtype> label_;
};

/**
* @brief Provides data to the Net from unified data file (text-image, image-image, text-text).
*/
template <typename Dtype>
class UnifiedDataLayer : public Layer<Dtype>, public InternalThread {
public:
    explicit UnifiedDataLayer(const LayerParameter& param);
    virtual ~UnifiedDataLayer() {}

    virtual inline const char* type() const { return "SiamesedData"; }
    virtual inline int ExactNumBottomBlobs() const { return 0; }

    virtual void Reshape(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top);
    virtual void DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void LayerSetUp(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top);
    void Forward_cpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top);
    void Forward_gpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
        Forward_cpu(bottom, top);
    }

    virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {}
    virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {}

    // Prefetches batches (asynchronously if to GPU memory)
    static const int PREFETCH_COUNT = 3;

    // Get attributes for libfeatureextractor
    UnifiedBatch<Dtype>* GetOneBatch() { return prefetch_full_.peek(); }
    vector<UnifiedDataLoader<Dtype>> data_loader_vec() { return data_loader_vec_; }

protected:
  // The thread's function
  virtual void InternalThreadEntry();
  virtual void load_batch(UnifiedBatch<Dtype>* batch);

  int batch_size_;
  size_t num_output_data_;

  UnifiedBatch<Dtype> prefetch_[PREFETCH_COUNT];
  BlockingQueue<UnifiedBatch<Dtype>*> prefetch_free_;
  BlockingQueue<UnifiedBatch<Dtype>*> prefetch_full_;
  vector<UnifiedDataLoader<Dtype>> data_loader_vec_;

  UnifiedDataParameter param_;
  vector<vector<string>> batch_samples_;
  vector<vector<ImageData>> batch_images_;

  UnifiedDataReader reader_;

  friend class UnifiedDataReader;
};

}  // namespace caffe

#endif  // CAFFE_UNIFIED_DATA_LAYER_HPP_
