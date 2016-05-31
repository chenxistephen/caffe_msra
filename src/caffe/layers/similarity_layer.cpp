#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>
#include <ppl.h>

#include "caffe/layers/similarity_layer.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe {

template <typename Dtype>
void SimilarityLayer<Dtype>::LayerSetUp(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
    Layer<Dtype>::LayerSetUp(bottom, top);
    auto& sb0 = bottom[0]->shape();
    auto& sb1 = bottom[1]->shape();
    CHECK((sb0.size() == 2 && sb1.size() == 2) || (sb0.size() == 4 && sb1.size() == 4 && sb0[2] == 1 && sb0[3] == 1 && sb1[2] == 1 && sb1[3] == 1))
        << "Input size incompatible for the two inputs.";
    CHECK(sb0[0] == sb1[0] && sb0[1] == sb1[1])
        << "Input size incompatible for the two inputs.";
}

template <typename Dtype>
void SimilarityLayer<Dtype>::Reshape(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
    auto& sb0 = bottom[0]->shape();
    auto& sb1 = bottom[1]->shape();
    CHECK((sb0.size() == 2 && sb1.size() == 2) || (sb0.size() == 4 && sb1.size() == 4 && sb0[2] == 1 && sb0[3] == 1 && sb1[2] == 1 && sb1[3] == 1))
        << "Input size incompatible for the two inputs.";
    CHECK(sb0[0] == sb1[0] && sb0[1] == sb1[1])
        << "Input size incompatible for the two inputs.";
    
    vector<int> top_shape{ sb0[0], 1 };
    top[0]->Reshape(top_shape);

}

template <typename Dtype>
void L2SimilarityLayer<Dtype>::Reshape(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
    SimilarityLayer<Dtype>::Reshape(bottom, top);

    vector<int> diff_shape{ bottom[0]->shape()[0], bottom[0]->shape()[1] };
    diff_.Reshape(diff_shape);
}

template <typename Dtype>
void L2SimilarityLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
    int count = bottom[0]->count();
    caffe_sub(
        count,
        bottom[0]->cpu_data(), 
        bottom[1]->cpu_data(),
        diff_.mutable_cpu_data());

    const int channels = bottom[0]->shape()[1];
    for (int i = 0; i < bottom[0]->num(); ++i) {
        top[0]->mutable_cpu_data()[i] = -1.0*caffe_cpu_dot(channels,
            diff_.cpu_data() + (i*channels), diff_.cpu_data() + (i*channels));
    }

}

template <typename Dtype>
void L2SimilarityLayer<Dtype>::Backward_cpu(
    const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    const Dtype* top_diff = top[0]->cpu_diff();
    for (int i = 0; i < 2; ++i) {
        if (propagate_down[i]) {
            Dtype* bottom_diff = bottom[i]->mutable_cpu_diff();
            const Dtype sign = (i == 0) ? -2. : 2.;

            int num = bottom[i]->shape()[0];
            int channels = bottom[i]->shape()[1];
            for (int j = 0; j < num; ++j) {
                const Dtype alpha = sign * top_diff[j];
                caffe_cpu_axpby(
                    channels,
                    alpha,
                    diff_.cpu_data() + (j*channels),
                    Dtype(0.0),
                    bottom_diff + (j*channels));
            }
        }
    }
}

template <typename Dtype>
void CosineSimilarityLayer<Dtype>::LayerSetUp(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
    SimilarityLayer<Dtype>::LayerSetUp(bottom, top);
    eps_ = this->layer_param_.cosine_similarity_param().eps();
}

template <typename Dtype>
void CosineSimilarityLayer<Dtype>::Reshape(
    const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
    SimilarityLayer<Dtype>::Reshape(bottom, top);

    vector<int> shape{ bottom[0]->shape()[0]};
    dotproduct_.Reshape(shape);
    norm0_.Reshape(shape);
    norm1_.Reshape(shape);
}


template <typename Dtype>
void CosineSimilarityLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
    const int num = bottom[0]->shape()[0];
    const int channels = bottom[0]->shape()[1];
    const Dtype* bottom0_data = bottom[0]->cpu_data();
    const Dtype* bottom1_data = bottom[1]->cpu_data();
    Dtype* top_data = top[0]->mutable_cpu_data();
    Dtype* dotproduct_data = dotproduct_.mutable_cpu_data();
    Dtype* norm0_data = norm0_.mutable_cpu_data();
    Dtype* norm1_data = norm1_.mutable_cpu_data();

    for (int i = 0; i < num; ++i) {
        dotproduct_data[i] = caffe_cpu_dot(channels,
            bottom0_data + (i*channels), bottom1_data + (i*channels));
        norm0_data[i] = caffe_cpu_dot(channels,
            bottom0_data + (i*channels), bottom0_data + (i*channels));
        norm1_data[i] = caffe_cpu_dot(channels,
            bottom1_data + (i*channels), bottom1_data + (i*channels));
    }

    // calculate cosine similarity
    caffe_powx(num, norm0_data, Dtype(0.5), norm0_data);
    caffe_powx(num, norm1_data, Dtype(0.5), norm1_data);

    caffe_mul(num, norm0_data, norm1_data, top_data);
    caffe_add_scalar(num, eps_, top_data);
    caffe_div(num, dotproduct_data, top_data, top_data);

    // scale it to [0,1]
    caffe_scal(num, Dtype(0.5), top_data);
    caffe_add_scalar(num, Dtype(0.5), top_data);

}

template <typename Dtype>
void CosineSimilarityLayer<Dtype>::Backward_cpu(
    const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    const Dtype* top_diff = top[0]->cpu_diff();
    const Dtype* bottom0_data = bottom[0]->cpu_data();
    const Dtype* bottom1_data = bottom[1]->cpu_data();
    const Dtype* dotproduct_data = dotproduct_.cpu_data();
    const Dtype* norm0_data = norm0_.cpu_data();
    const Dtype* norm1_data = norm1_.cpu_data();
    Dtype* bottom0_diff = bottom[0]->mutable_cpu_diff();
    Dtype* bottom1_diff = bottom[1]->mutable_cpu_diff();
    
    const int num = bottom[0]->shape()[0];
    const int channels = bottom[0]->shape()[1];

    for (int i = 0; i < num; ++i) {	

        if (propagate_down[0]) {
            Dtype alpha = 0.5 / std::max(norm0_data[i] * norm1_data[i], eps_);
            Dtype beta = -0.5 * dotproduct_data[i] / std::max(pow(norm0_data[i], Dtype(3.0)) * norm1_data[i], eps_);

            caffe_copy(channels, bottom0_data + (i*channels), bottom0_diff + (i*channels));
            caffe_cpu_axpby(channels, alpha, bottom1_data + (i*channels), beta, bottom0_diff + (i*channels));
            caffe_scal(channels, top_diff[i], bottom0_diff + (i*channels));
        }

        if (propagate_down[1]) {
            Dtype alpha = 0.5 / std::max(norm0_data[i] * norm1_data[i], eps_);
            Dtype beta = -0.5 * dotproduct_data[i] / std::max(pow(norm1_data[i], Dtype(3.0)) * norm0_data[i], eps_);

            caffe_copy(channels, bottom1_data + (i*channels), bottom1_diff + (i*channels));
            caffe_cpu_axpby(channels, alpha, bottom0_data + (i*channels), beta, bottom1_diff + (i*channels));
            caffe_scal(channels, top_diff[i], bottom1_diff + (i*channels));
        }
    }

}


#ifdef CPU_ONLY
STUB_GPU(L2SimilarityLayer);
STUB_GPU(CosineSimilarityLayer);
#endif

INSTANTIATE_CLASS(L2SimilarityLayer);
REGISTER_LAYER_CLASS(L2Similarity);
INSTANTIATE_CLASS(CosineSimilarityLayer);
REGISTER_LAYER_CLASS(CosineSimilarity);

}  // namespace caffe
