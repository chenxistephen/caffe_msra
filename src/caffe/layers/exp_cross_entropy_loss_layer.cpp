#include <algorithm>
#include <cfloat>
#include <vector>

#include "caffe/layer.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/layers/exp_cross_entropy_loss_layer.hpp"

namespace caffe {

template <typename Dtype>
void ExpCrossEntropyLossLayer<Dtype>::LayerSetUp(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
    LossLayer<Dtype>::LayerSetUp(bottom, top);
    eps_ = this->layer_param_.exp_cross_entropy_loss_param().eps();
}

template <typename Dtype>
void ExpCrossEntropyLossLayer<Dtype>::Reshape(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  LossLayer<Dtype>::Reshape(bottom, top);
  CHECK_EQ(bottom[0]->count(), bottom[1]->count()) <<
      "EXP_CROSS_ENTROPY_LOSS layer inputs must have the same count.";
  exp_output_->Reshape(bottom[0]->shape());
}

template <typename Dtype>
void ExpCrossEntropyLossLayer<Dtype>::Forward_cpu(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  // Compute the loss (negative log likelihood)
  const int count = bottom[0]->count();
  const int num = bottom[0]->num();
  const Dtype* input_data = bottom[0]->cpu_data();
  const Dtype* target = bottom[1]->cpu_data();
  Dtype* exp_output = exp_output_->mutable_cpu_data();
  Dtype loss = 0;
  caffe_exp(count, input_data, exp_output);
  for (int i = 0; i < count; ++i) {
      loss -= target[i] * input_data[i] + (1.0 - target[i]) * log(std::max(Dtype(1.0) - exp_output[i], eps_));
  }
  top[0]->mutable_cpu_data()[0] = loss / num;
}

template <typename Dtype>
void ExpCrossEntropyLossLayer<Dtype>::Backward_cpu(
    const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
  if (propagate_down[1]) {
    LOG(FATAL) << this->type()
               << " Layer cannot backpropagate to label inputs.";
  }
  if (propagate_down[0]) {
    // First, compute the diff
    const int count = bottom[0]->count();
    const int num = bottom[0]->num();
    const Dtype* input_data = bottom[0]->cpu_data();
    const Dtype* target = bottom[1]->cpu_data();
    const Dtype* exp_output = exp_output_->cpu_data();
    Dtype* bottom_diff = bottom[0]->mutable_cpu_diff();
    for (int i = 0; i < count; ++i) {
        bottom_diff[i] = -target[i] + (1.0 - target[i]) / std::max(Dtype(1.0) - exp_output[i], eps_) * exp_output[i];
    }
    // Scale down gradient
    const Dtype loss_weight = top[0]->cpu_diff()[0];
    caffe_scal(count, loss_weight / num, bottom_diff);
  }
}

#ifdef CPU_ONLY
STUB_GPU(ExpCrossEntropyLossLayer);
#endif

INSTANTIATE_CLASS(ExpCrossEntropyLossLayer);
REGISTER_LAYER_CLASS(ExpCrossEntropyLoss);

}  // namespace caffe
