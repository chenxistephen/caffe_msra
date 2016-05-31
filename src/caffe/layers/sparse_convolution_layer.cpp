#include "caffe/layers/nlp_layers.hpp"
#include "caffe/filler.hpp"
#include <ppl.h>
using namespace concurrency;

namespace caffe{

template<typename Dtype>
void SparseConvolutionLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top){
    CHECK(bottom.size() >= 1 && bottom[0]->SparseBlob().get() != nullptr)
        << "The input should be sparse blob"; 

    SparseConvolutionParameter sparse_conv_param = this->layer_param_.sparse_conv_param();
    kernel_size_ = sparse_conv_param.kernel_size();
    stride_ = sparse_conv_param.stride();
    num_outputs_ = sparse_conv_param.num_outputs();
    bias_term_ = sparse_conv_param.bias_term();
    vocab_size_ = bottom[0]->SparseBlob()->vocab_size();
    share_weights_inside_kernel_ = sparse_conv_param.share_weights_inside_kernel();
    weight_dim_ = share_weights_inside_kernel_ ? (int)vocab_size_ : (int)vocab_size_ * kernel_size_;

    if (this->blobs_.size() > 0) {
        LOG(INFO) << "Skipping parameter initialization";
    }
    else {
        if (bias_term_) {
            this->blobs_.resize(2);
        }
        else {
            this->blobs_.resize(1);
        }
        // Initialize and fill the weights:
        vector<int> weight_shape{ (int)num_outputs_, weight_dim_ };
        this->blobs_[0].reset(new Blob<Dtype>(weight_shape));
        shared_ptr<Filler<Dtype> > weight_filler(GetFiller<Dtype>(
            this->layer_param_.sparse_conv_param().weight_filler()));
        weight_filler->Fill(this->blobs_[0].get());
        // If necessary, initialize and fill the biases.
        if (bias_term_) {
            vector<int> bias_shape{ 1, (int)num_outputs_ };
            this->blobs_[1].reset(new Blob<Dtype>(bias_shape));
            shared_ptr<Filler<Dtype> > bias_filler(GetFiller<Dtype>(
                this->layer_param_.sparse_conv_param().bias_filler()));
            bias_filler->Fill(this->blobs_[1].get());
        }

        // Propagate gradients to the parameters (as directed by backward pass).
        this->param_propagate_down_.resize(this->blobs_.size(), true);
    }
}

template<typename Dtype>
void SparseConvolutionLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top){
    CHECK(top.size()>=1)
        << "There is at least one output from SparseConvolutionLayer.";
    CHECK(bottom.size()>=1)
        << "There is at least one input from SparseConvolutionLayer.";
    int count_images = bottom[0]->SparseBlob()->batch_size();
    auto word_count = bottom[0]->SparseBlob()->word_count();
    int out_height = (word_count - 1) / stride_ + 1;
    for (int i = 0; i < bottom.size(); ++i) {
        vector<int> top_shape(4);
        top_shape[0] = count_images;
        top_shape[1] = num_outputs_;
        top_shape[2] = out_height;
        top_shape[3] = 1;
        top[i]->Reshape(top_shape);
    }
    if (bias_term_) {
        bias_multiplier_a_.Reshape(vector<int> { 1, out_height });
        caffe_set(bias_multiplier_a_.count(), Dtype(1.0), bias_multiplier_a_.mutable_cpu_data());
        bias_tmp_.Reshape(vector<int> { (int)num_outputs_, out_height });
        bias_multiplier_b_.Reshape(vector<int> { 1, count_images });
        caffe_set(bias_multiplier_b_.count(), Dtype(1), bias_multiplier_b_.mutable_cpu_data());
    }
}

template<typename Dtype>
void SparseConvolutionLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
    const vector<Blob<Dtype>*>& top){
    for (int i = 0; i < bottom.size(); ++i) {
        CHECK(bottom[i]->SparseBlob().get() != nullptr)
            << "The input should be sparse blob";
        auto sparse_input = bottom[i]->SparseBlob();
        auto& weights = *this->blobs_[0];

        auto num_samples = sparse_input->batch_size();
        auto word_count = sparse_input->word_count();
        auto out_height = (word_count - 1) / stride_ + 1;

        //Blob<Dtype> tmp_output(vector < int > {(int)num_samples, (int)num_outputs_, (int)word_count});

        caffe_set(top[i]->count(), Dtype(0.0), top[i]->mutable_cpu_data());
        const SparseItem<Dtype>* p_sprase_items = sparse_input->cpu_data();
        int valid_count = GetNonZeroCount(*sparse_input);
        auto top_data = top[i];
        Dtype* top_cpu_data = top_data->mutable_cpu_data();
        const Dtype* weights_cpu = weights.cpu_data();
        for (int sid = 0; sid < valid_count; sid++)
        {
            //if (!p_sprase_items[sid].IsNull()){
            for (int out_id = p_sprase_items[sid].word_id / stride_; out_id >= 0; --out_id){
                int offset = p_sprase_items[sid].word_id - out_id*stride_;
                if (offset < kernel_size_){
                    const Dtype* weights_cpu_offset = share_weights_inside_kernel_ ? weights_cpu: (weights_cpu + offset * vocab_size_);
                    Dtype* p_out = top_cpu_data + top_data->offset(p_sprase_items[sid].sample_id, 0, out_id);

                    for (auto k = 0; k < num_outputs_; ++k){
                        Dtype new_val = p_sprase_items[sid].val * weights_cpu_offset[k*weight_dim_ + p_sprase_items[sid].feature_id];
                        p_out[k*out_height] += new_val;
                    }
                }
                else
                    break;
            }
        }

        if (this->bias_term_) {
            const Dtype* bias = this->blobs_[1]->cpu_data();
            caffe_cpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, num_outputs_,
                out_height, 1, (Dtype)1., bias, bias_multiplier_a_.cpu_data(),
                (Dtype)0., bias_tmp_.mutable_cpu_data());
            caffe_cpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, num_samples,
                num_outputs_*out_height, 1, (Dtype)1., bias_multiplier_b_.cpu_data(), bias_tmp_.cpu_data(),
                (Dtype)1., top[i]->mutable_cpu_data());
        }
    }
}

template<class Dtype>
int SparseConvolutionLayer<Dtype>::GetNonZeroCount(const SparseBlob<Dtype>& sparse_input){
    int count = 0;
    const SparseItem<Dtype>* p_sprase_items = sparse_input.cpu_data();
    for (int i = 0; i < sparse_input.count(); ++i, ++count){
        if (p_sprase_items[i].IsNull()) break;
    }
    return count;
}


template<typename Dtype>
void SparseConvolutionLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom){
    const Dtype* weight = this->blobs_[0]->cpu_data();
    Dtype* weight_diff = this->blobs_[0]->mutable_cpu_diff();
    if (this->param_propagate_down_[0]) {
        caffe_set(this->blobs_[0]->count(), Dtype(0), weight_diff);
    }
    if (this->bias_term_ && this->param_propagate_down_[1]) {
        caffe_set(this->blobs_[1]->count(), Dtype(0),
            this->blobs_[1]->mutable_cpu_diff());
    }
    for (int i = 0; i < top.size(); ++i) {
        CHECK(bottom[i]->SparseBlob().get() != nullptr)
            << "The input should be sparse blob";
        auto top_data = top[i];
        const Dtype* top_diff = top_data->cpu_diff();
        auto& sparse_input = bottom[i]->SparseBlob();
        

        auto& weights = *this->blobs_[0];
        auto num_samples = sparse_input->batch_size();
        auto word_count = sparse_input->word_count();
        auto out_height = (word_count - 1) / stride_ + 1;
        const SparseItem<Dtype>* p_sprase_items = sparse_input->cpu_data();

        // Bias gradient, if necessary.
        if (this->bias_term_ && this->param_propagate_down_[1]) {
            Dtype* bias_diff = this->blobs_[1]->mutable_cpu_diff();
            caffe_cpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, 1, 
                num_outputs_* out_height, num_samples, (Dtype)1., bias_multiplier_b_.cpu_data(), top_diff,
                (Dtype)0., bias_tmp_.mutable_cpu_data());

            caffe_cpu_gemv<Dtype>(CblasNoTrans, num_outputs_, out_height, 1.,
                bias_tmp_.mutable_cpu_data(), bias_multiplier_a_.cpu_data(), 0., bias_diff);
        }
        if (this->param_propagate_down_[0]) {
            int valid_count = GetNonZeroCount(*sparse_input);
            for (int sid = 0; sid < valid_count; sid++) {
                int nSample = p_sprase_items[sid].sample_id;
                for (int out_id = p_sprase_items[sid].word_id / stride_; out_id >= 0; --out_id){
                    const Dtype* p_in = top_diff + top_data->offset(p_sprase_items[sid].sample_id, 0, out_id);
                    int offset = p_sprase_items[sid].word_id - out_id*stride_;
                    if (offset < kernel_size_){
                        Dtype* weights_diff_offset = share_weights_inside_kernel_ ? weight_diff:  (weight_diff + offset * vocab_size_);
                        for (auto k = 0; k < num_outputs_; ++k){
                            Dtype new_val = p_sprase_items[sid].val * p_in[k*out_height];
                            weights_diff_offset[k*weight_dim_ + p_sprase_items[sid].feature_id] += new_val;
                        }
                    }
                }
            }
        }
        CHECK(!propagate_down[i]) << "there is no propagation down for sparse convolution layer";
    }
}

#ifdef CPU_ONLY
STUB_GPU(SparseConvolutionLayer);
#endif

INSTANTIATE_CLASS(SparseConvolutionLayer);
REGISTER_LAYER_CLASS(SparseConvolution);
}