#ifndef CAFFE_EXP_CROSS_ENTROPY_LOSS_LAYER_HPP_
#define CAFFE_EXP_CROSS_ENTROPY_LOSS_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "caffe/layers/loss_layer.hpp"

namespace caffe {

template <typename Dtype>
class ExpCrossEntropyLossLayer : public LossLayer<Dtype> {
public:
  explicit ExpCrossEntropyLossLayer(const LayerParameter& param)
    : LossLayer<Dtype>(param),
    exp_output_(new Blob<Dtype>()) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "ExpCrossEntropyLoss"; }

protected:
  /// @copydoc ExpCrossEntropyLossLayer
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);

  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

  shared_ptr<Blob<Dtype> > exp_output_;
  Dtype eps_;

};

}  // namespace caffe

#endif  // CAFFE_EXP_CROSS_ENTROPY_LOSS_LAYER_HPP_
