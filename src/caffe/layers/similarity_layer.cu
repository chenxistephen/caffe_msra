#include "caffe/layers/similarity_layer.hpp"

namespace caffe {

template <typename Dtype>
void L2SimilarityLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
    Forward_cpu(bottom, top);
}

template <typename Dtype>
void L2SimilarityLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    Backward_cpu(top, propagate_down, bottom);
}

template <typename Dtype>
void CosineSimilarityLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
    Forward_cpu(bottom, top);
}

template <typename Dtype>
void CosineSimilarityLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    Backward_cpu(top, propagate_down, bottom);
}

INSTANTIATE_LAYER_GPU_FUNCS(L2SimilarityLayer);
INSTANTIATE_LAYER_GPU_FUNCS(CosineSimilarityLayer);

};