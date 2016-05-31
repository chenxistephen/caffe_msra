#include <algorithm>
#include <cfloat>
#include <vector>

#include "caffe/common.hpp"
#include "caffe/layer.hpp"
#include "caffe/syncedmem.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/layers/affine_transform_layer.hpp"

namespace caffe {

using std::min;
using std::max;

template <typename Dtype>
void AffineMatConversionLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
}

template <typename Dtype>
void AffineMatConversionLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
    top[0]->Reshape(vector < int > {bottom[0]->shape()[0], 6});
    const auto& input_dim = bottom[0]->shape();
    int input_mat_len = input_dim[1];
    CHECK_EQ(bottom[0]->count(2), 1) << "input dimension is incorrect";

    const auto& conversion_param = this->layer_param_.affine_mat_conversion_param();
    switch (conversion_param.conversion_type()) {
    case AffineMatConversionParameter::RECT_CROP:
        CHECK_EQ(input_mat_len, 4) << "RECT_CROP should have 4-dim transformation mat";
        break;
    default:
        LOG(FATAL) << "incorrect conversion type";
    }
}


template <typename Dtype>
void AffineMatConversionLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
    const auto& conversion_param = this->layer_param_.affine_mat_conversion_param();

    const auto& input_dim = bottom[0]->shape();
    int input_mat_len = input_dim[1];

    for (int i = 0; i < input_dim[0]; ++i) {
        const Dtype* p_in = bottom[0]->cpu_data() + input_mat_len*i;
        Dtype* p_out = top[0]->mutable_cpu_data() + 6*i;
        switch (conversion_param.conversion_type()) {
        case AffineMatConversionParameter::RECT_CROP:
            p_out[0] = p_in[0]; p_out[1] = 0;       p_out[2] = p_in[1];
            p_out[3] = 0;       p_out[4] = p_in[2]; p_out[5] = p_in[3];
            break;
        default:
            LOG(FATAL) << "incorrect conversion type";
        }
    }
}

template <typename Dtype>
void AffineMatConversionLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    if (!propagate_down[0]) {
        return;
    }
    const auto& conversion_param = this->layer_param_.affine_mat_conversion_param();

    const auto& input_dim = bottom[0]->shape();
    int input_mat_len = input_dim[1];

    for (int i = 0; i < input_dim[0]; ++i) {
        const Dtype* p_top_diff = top[0]->cpu_diff() + 6*i;
        Dtype* p_bottom_diff = bottom[0]->mutable_cpu_diff() + input_mat_len * i;
        switch (conversion_param.conversion_type()) {
        case AffineMatConversionParameter::RECT_CROP:
            p_bottom_diff[0] = p_top_diff[0]; p_bottom_diff[1] = p_top_diff[2]; p_bottom_diff[2] = p_top_diff[4]; p_bottom_diff[3] = p_top_diff[5];
            break;
        default:
            LOG(FATAL) << "incorrect conversion type";
        }
    }
}
template <typename Dtype>
void AffineMatConversionLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top) {
    Forward_cpu(bottom, top);
}

template <typename Dtype>
void AffineMatConversionLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    Backward_cpu(top, propagate_down, bottom);
}

/*
TODO: current gpu code is run in cpu, we need to implement the gpu part and uncomment here.
#ifdef CPU_ONLY
STUB_GPU(AffineMatConversionLayer);
#endif
*/

INSTANTIATE_CLASS(AffineMatConversionLayer);
REGISTER_LAYER_CLASS(AffineMatConversion);

}  // namespace caffe
