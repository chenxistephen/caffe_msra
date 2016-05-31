#include <algorithm>
#include <cfloat>
#include <vector>
#include "caffe/layer.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/layers/nlp_layers.hpp"

namespace caffe {
    template<typename Dtype>
    __global__ void SparseConvKernel(const SparseItem<Dtype>* p_sprase_items, const Dtype* weights_gpu, Dtype* top_gpu_data, const int valid_count, const int vocab_size, const int kernel_size, const int stride, const bool share_weights_inside_kernel, const int weight_dim, const int out_height, const int out_feadim);

    template<>
    __global__ void SparseConvKernel<float>(const SparseItem<float>* p_sprase_items, const float* weights_gpu, float* top_gpu_data, const int valid_count, const int vocab_size, const int kernel_size, const int stride, const bool share_weights_inside_kernel, const int weight_dim, const int out_height, const int out_feadim) {
        uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
        uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
        if (idx < valid_count && idy < out_feadim) {
            for (int out_id = p_sprase_items[idx].word_id / stride; out_id >= 0; --out_id){
                int offset = p_sprase_items[idx].word_id - out_id*stride;
                if (offset < kernel_size){
                    const float* weights_gpu_offset = share_weights_inside_kernel ? weights_gpu : (weights_gpu + offset * vocab_size);
                    float* p_out = top_gpu_data + p_sprase_items[idx].sample_id * out_feadim * out_height + out_id;
                    float new_val = p_sprase_items[idx].val * weights_gpu_offset[idy * weight_dim + p_sprase_items[idx].feature_id];
                    atomicAdd(p_out + idy * out_height, new_val);
                }
                else {
                    break;
                }
            }
        }
    }
    
    __device__ inline void atomicAdd_double(double *address, double value)
    {
        unsigned long long oldval, newval, readback;

        oldval = __double_as_longlong(*address);
        newval = __double_as_longlong(__longlong_as_double(oldval) + value);
        while ((readback = atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
        {
            oldval = readback;
            newval = __double_as_longlong(__longlong_as_double(oldval) + value);
        }
    }

    template<>
    __global__ void SparseConvKernel<double>(const SparseItem<double>* p_sprase_items, const double* weights_gpu, double* top_gpu_data, const int valid_count, const int vocab_size, const int kernel_size, const int stride, const bool share_weights_inside_kernel, const int weight_dim, const int out_height, const int out_feadim) {
        uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
        uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
        if (idx < valid_count && idy < out_feadim) {
            for (int out_id = p_sprase_items[idx].word_id / stride; out_id >= 0; --out_id){
                int offset = p_sprase_items[idx].word_id - out_id*stride;
                if (offset < kernel_size){
                    const double* weights_gpu_offset = share_weights_inside_kernel ? weights_gpu : (weights_gpu + offset * vocab_size);
                    double* p_out = top_gpu_data + p_sprase_items[idx].sample_id * out_feadim * out_height + out_id;
                    double new_val = p_sprase_items[idx].val * weights_gpu_offset[idy * weight_dim + p_sprase_items[idx].feature_id];
                    atomicAdd_double(p_out + idy * out_height, new_val);
                }
                else {
                    break;
                }
            }
        }
    }

    template <typename Dtype>
    void SparseConvolutionLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
        for (int i = 0; i < bottom.size(); ++i) {
            CHECK(bottom[i]->SparseBlob().get() != nullptr)
                << "The input should be sparse blob";
            auto sparse_input = bottom[i]->SparseBlob();
            int valid_count = GetNonZeroCount(*sparse_input);
            auto& weights = *this->blobs_[0];

            auto num_samples = sparse_input->batch_size();
            auto word_count = sparse_input->word_count();
            auto out_height = (word_count - 1) / stride_ + 1;

            caffe_gpu_set(top[i]->count(), Dtype(0.0), top[i]->mutable_gpu_data());
            const SparseItem<Dtype>* p_sparse_items = sparse_input->gpu_data();
            auto top_data = top[i];
            Dtype* top_gpu_data = top_data->mutable_gpu_data();
            const Dtype* weights_gpu = weights.gpu_data();

            dim3 thread_tail(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
            dim3 block_tail((valid_count + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (num_outputs_ + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
            SparseConvKernel<Dtype> << <block_tail, thread_tail >> >(p_sparse_items, weights_gpu, top_gpu_data, valid_count, vocab_size_, kernel_size_, stride_, share_weights_inside_kernel_, weight_dim_, out_height, num_outputs_);
            //SparseConvKernel<Dtype> << <CAFFE_GET_BLOCKS(valid_count), CAFFE_CUDA_NUM_THREADS >> >(p_sparse_items, weights_gpu, top_gpu_data, valid_count, vocab_size_, kernel_size_, stride_, share_weights_inside_kernel_, weight_dim_, out_height, num_outputs_);

            if (this->bias_term_) {
                const Dtype* bias = this->blobs_[1]->gpu_data();
                caffe_gpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, num_outputs_,
                    out_height, 1, (Dtype)1., bias, bias_multiplier_a_.gpu_data(),
                    (Dtype)0., bias_tmp_.mutable_gpu_data());

                caffe_gpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, num_samples,
                    num_outputs_*out_height, 1, (Dtype)1., bias_multiplier_b_.gpu_data(), bias_tmp_.gpu_data(),
                    (Dtype)1., top[i]->mutable_gpu_data());
            }
            
        }
    }

    template<typename Dtype>
    __global__ void SparseConvBackKernel(const SparseItem<Dtype>* p_sprase_items, const Dtype* top_diff, Dtype* weight_diff, const int valid_count, const int vocab_size, const int kernel_size, const int stride, const bool share_weights_inside_kernel, const int weight_dim, const int out_height, const int out_feadim);

    template<>
    __global__ void SparseConvBackKernel<float>(const SparseItem<float>* p_sprase_items, const float* top_diff, float* weight_diff, const int valid_count, const int vocab_size, const int kernel_size, const int stride, const bool share_weights_inside_kernel, const int weight_dim, const int out_height, const int out_feadim) {
        uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
        uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
        if (idx < valid_count && idy < out_feadim) {
            for (int out_id = p_sprase_items[idx].word_id / stride; out_id >= 0; --out_id){
                const float* p_in = top_diff + p_sprase_items[idx].sample_id * out_feadim * out_height + out_id;
                int offset = p_sprase_items[idx].word_id - out_id * stride;
                if (offset < kernel_size){
                    float* weights_diff_offset = share_weights_inside_kernel ? weight_diff : (weight_diff + offset * vocab_size);
                    float new_val = p_sprase_items[idx].val * p_in[idy*out_height];
                    atomicAdd(weights_diff_offset + (idy * weight_dim + p_sprase_items[idx].feature_id), new_val);
                }
                else {
                    break;
                }
            }
        }
    }


    template<>
    __global__ void SparseConvBackKernel<double>(const SparseItem<double>* p_sprase_items, const double* top_diff, double* weight_diff, const int valid_count, const int vocab_size, const int kernel_size, const int stride, const bool share_weights_inside_kernel, const int weight_dim, const int out_height, const int out_feadim) {
        uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
        uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
        if (idx < valid_count && idy < out_feadim) {
            for (int out_id = p_sprase_items[idx].word_id / stride; out_id >= 0; --out_id){
                const double* p_in = top_diff + p_sprase_items[idx].sample_id * out_feadim * out_height + out_id;
                int offset = p_sprase_items[idx].word_id - out_id * stride;
                if (offset < kernel_size){
                    double* weights_diff_offset = share_weights_inside_kernel ? weight_diff : (weight_diff + offset * vocab_size);
                    double new_val = p_sprase_items[idx].val * p_in[idy*out_height];
                    atomicAdd_double(weights_diff_offset + (idy * weight_dim + p_sprase_items[idx].feature_id), new_val);
                }
                else {
                    break;
                }
            }
        }
    }

    template<typename Dtype>
    void SparseConvolutionLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom){
        const Dtype* weight = this->blobs_[0]->gpu_data();
        Dtype* weight_diff = this->blobs_[0]->mutable_gpu_diff();
        if (this->param_propagate_down_[0]) {
            caffe_gpu_set(this->blobs_[0]->count(), Dtype(0), weight_diff);
        }
        if (this->bias_term_ && this->param_propagate_down_[1]) {
            caffe_gpu_set(this->blobs_[1]->count(), Dtype(0),
                this->blobs_[1]->mutable_gpu_diff());
        }
        for (int i = 0; i < top.size(); ++i) {
            CHECK(bottom[i]->SparseBlob().get() != nullptr)
                << "The input should be sparse blob";
            auto top_data = top[i];
            const Dtype* top_diff = top_data->gpu_diff();
            auto& sparse_input = bottom[i]->SparseBlob();
            int valid_count = GetNonZeroCount(*sparse_input);

            auto& weights = *this->blobs_[0];
            auto num_samples = sparse_input->batch_size();
            auto word_count = sparse_input->word_count();
            auto out_height = (word_count - 1) / stride_ + 1;
            const SparseItem<Dtype>* p_sprase_items = sparse_input->gpu_data();

            
            // Bias gradient, if necessary.
            if (this->bias_term_ && this->param_propagate_down_[1]) {
                Dtype* bias_diff = this->blobs_[1]->mutable_gpu_diff();
                caffe_gpu_gemm<Dtype>(CblasNoTrans, CblasNoTrans, 1,
                    num_outputs_* out_height, num_samples, (Dtype)1., bias_multiplier_b_.gpu_data(), top_diff,
                    (Dtype)0., bias_tmp_.mutable_gpu_data());

                caffe_gpu_gemv<Dtype>(CblasNoTrans, num_outputs_, out_height, 1.,
                    bias_tmp_.mutable_gpu_data(), bias_multiplier_a_.gpu_data(), 0., bias_diff);
            }
            
            if (this->param_propagate_down_[0]) {
                dim3 thread_tail(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
                dim3 block_tail((valid_count + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (num_outputs_ + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
                SparseConvBackKernel<Dtype> << <block_tail, thread_tail >> >(p_sprase_items, top_diff, weight_diff, valid_count, vocab_size_, kernel_size_, stride_, share_weights_inside_kernel_, weight_dim_, out_height, num_outputs_);
            }
            CHECK(!propagate_down[i], "there is no propagation down for sparse convolution layer");
            
        }
    }

    INSTANTIATE_LAYER_GPU_FUNCS(SparseConvolutionLayer);
}