#ifndef CAFFE_NLP_LAYERS_HPP_
#define CAFFE_NLP_LAYERS_HPP_

#include <string>
#include <utility>
#include <vector>

#include "caffe/blob.hpp"
#include "caffe/common.hpp"
#include "caffe/layer.hpp"
#include "caffe/proto/caffe.pb.h"

namespace caffe {

template <typename Dtype>
class SparseConvolutionLayer : public Layer < Dtype > {
public:
    explicit SparseConvolutionLayer(const LayerParameter& param)
        : Layer<Dtype>(param) {}

    void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void Reshape(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);

protected:
    void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
    void Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
    int GetNonZeroCount(const SparseBlob<Dtype>& sparse_input);

protected:
    uint32_t kernel_size_;
    uint32_t stride_;
    uint32_t num_outputs_;
    int vocab_size_;
    bool bias_term_;
    Blob<Dtype> bias_multiplier_a_;
    Blob<Dtype> bias_multiplier_b_;
    Blob<Dtype> bias_tmp_;
    bool share_weights_inside_kernel_;
    int weight_dim_;
};


template <typename Dtype>
class Sparse2DenseBlobConversionLayer : public Layer < Dtype > {
public:
    explicit Sparse2DenseBlobConversionLayer(const LayerParameter& param)
        : Layer<Dtype>(param) {}

    void LayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void Reshape(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);

protected:
    void Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top);
    void Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);
    void Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom);

protected:
    uint32_t kernel_size_;
    uint32_t stride_;
    uint32_t num_outputs_;
    bool bias_term_;
    Blob<Dtype> bias_multiplier_;
};

}

#endif
