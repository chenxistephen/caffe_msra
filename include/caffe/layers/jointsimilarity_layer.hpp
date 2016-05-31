#ifndef CAFFE_JOINTSIMILARITY_LAYER_HPP_
#define CAFFE_JOINTSIMILARITY_LAYER_HPP_

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
class JointSimilarityLayer : public Layer<Dtype> {
public:
  explicit JointSimilarityLayer(const LayerParameter& param)
    : Layer<Dtype>(param) {}
  virtual void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);

  virtual inline int ExactNumBottomBlobs() const { return 2; }
  virtual inline int ExactNumTopBlobs() const { return 1; }

protected:
  virtual void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  virtual Dtype ComputeSimilarity_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, int dim, const Dtype eps) = 0;
  virtual void PropagateDown_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, Dtype top_diff, Dtype* bottom_diff_row, int dim, const Dtype eps) = 0;
  void Sample_pairs(const int num, const int dim);

protected:
  shared_ptr<Caffe::RNG> column_sample_rng_;
  Dtype eps_;
  int column_sample_num_;
  // sample_mat0_[i, j] means the (sample_mat0_[i, j]-1) th row of bottom1 has been sampled out to compute similarity with ith row of bottom0
  // for j >= column_sample_num, sample_mat0_[i, j] should be equal to 0
  // size: num * num
  Blob<Dtype> sample_mat0_;
  // sample_mat1_[i, j] means the ith row of bottom1 has been sampled out to compute the (sample_mat1_[i, j] - 1) th value of top matrix, which means, it use the (sample_mat1_[i, j] - 1) / column_sample_num_ th row of bottom0 is used
  // if sample_mat1_[i, j] == 0, it means the ith row of bottom1 has been computed for less than j times, it is used for control for loop
  // size: num * num;
  Blob<Dtype> sample_mat1_;
  // the space to keep error derivative for each pair before summing them up
  Blob<Dtype> derivative_mat0_, derivative_mat1_;
  bool pairwise_;
};

/**
* @brief Compute the L2 similarity of the inputs.
*
* TODO(dox): thorough documentation for Forward, Backward, and proto params.
*/
template <typename Dtype>
class L2JointSimilarityLayer : public JointSimilarityLayer<Dtype> {
public:
  explicit L2JointSimilarityLayer(const LayerParameter& param)
    : JointSimilarityLayer<Dtype>(param) {}
  virtual inline const char* type() const { return "L2JointSimilarity"; }

protected:
  virtual Dtype ComputeSimilarity_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, int dim, const Dtype eps);
  virtual void PropagateDown_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, Dtype top_diff, Dtype* bottom_diff_row, int dim, const Dtype eps);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
};

/**
* @brief Compute the cos similarity of the inputs.
*
* TODO(dox): thorough documentation for Forward, Backward, and proto params.
*/
template <typename Dtype>
class CosineJointSimilarityLayer : public JointSimilarityLayer<Dtype> {
public:
  explicit CosineJointSimilarityLayer(const LayerParameter& param)
    : JointSimilarityLayer<Dtype>(param) {}
  virtual inline const char* type() const { return "CosineSimilarity"; }
  //virtual void Reshape(const vector<Blob<Dtype>*>& bottom,
  //    const vector<Blob<Dtype>*>& top);

protected:
  virtual Dtype ComputeSimilarity_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, int dim, const Dtype eps);
  virtual void PropagateDown_cpu(const Dtype* bottom_row0, const Dtype* bottom_row1, Dtype top_diff, Dtype* bottom_diff_row, int dim, const Dtype eps);
  virtual void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top);
  virtual void Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
  Blob<Dtype> norm0_;
  Blob<Dtype> norm1_;
  Blob<Dtype> dotproduct_;
};

}  // namespace caffe

#endif  // CAFFE_JOINTSIMILARITY_LAYER_HPP_
