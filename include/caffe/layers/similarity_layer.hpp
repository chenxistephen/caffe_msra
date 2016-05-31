#ifndef CAFFE_SIMILARITY_LAYER_HPP_
#define CAFFE_SIMILARITY_LAYER_HPP_

#include <vector>

#include "caffe/blob.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

#include "caffe/layers/neuron_layer.hpp"

namespace caffe {

/**
* @brief Compute the similarity of the inputs.
*
* TODO(dox): thorough documentation for Forward, Backward, and proto params.
*/
template <typename Dtype>
class SimilarityLayer : public Layer<Dtype> {
public:
  explicit SimilarityLayer(const LayerParameter& param)
    : Layer<Dtype>(param) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);

  virtual inline int ExactNumBottomBlobs() const { return 2; }
  virtual inline int ExactNumTopBlobs() const { return 1; }
};

/**
* @brief Compute the L2 similarity of the inputs.
*
* TODO(dox): thorough documentation for Forward, Backward, and proto params.
*/
template <typename Dtype>
class L2SimilarityLayer : public SimilarityLayer<Dtype> {
public:
  explicit L2SimilarityLayer(const LayerParameter& param)
    : SimilarityLayer<Dtype>(param) {}
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "L2Similarity"; }

protected:
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

private:
  Blob<Dtype> diff_;
};

/**
* @brief Compute the cos similarity of the inputs.
*
* TODO(dox): thorough documentation for Forward, Backward, and proto params.
*/
template <typename Dtype>
class CosineSimilarityLayer : public SimilarityLayer<Dtype> {
public:
  explicit CosineSimilarityLayer(const LayerParameter& param)
    : SimilarityLayer<Dtype>(param) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);

  virtual inline const char* type() const { return "CosineSimilarity"; }

protected:
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

private:
  Blob<Dtype> dotproduct_;
  Blob<Dtype> norm0_;
  Blob<Dtype> norm1_;
  Dtype eps_;
};

}  // namespace caffe

#endif  // CAFFE_SIMILARITY_LAYER_HPP_
