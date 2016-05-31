#include "caffe/layers/jointsimilarity_layer.hpp"

namespace caffe {
template <typename Dtype>
__global__ void L2DistKernel(const int num, const int dim, const int column_sample_num, const Dtype* bottom_data0, const Dtype* bottom_data1, const Dtype* sample_mat, Dtype* dist) {
    uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num && idy < column_sample_num) {
        const Dtype* bottom_row0 = bottom_data0 + idx * dim;
        const Dtype* bottom_row1 = bottom_data1 + int(sample_mat[idx * num + idy] - 1) * dim;
        Dtype sum = 0;
        for (int j = 0; j < dim; ++j) {
            Dtype d = bottom_row0[j] - bottom_row1[j];
            sum -= d*d;
        }
        dist[idx * column_sample_num + idy] = sum;
    }
}

template <typename Dtype>
__global__ void CosineDistKernel(const int num, const int dim, const int column_sample_num, const Dtype* bottom_data0, const Dtype* bottom_data1, const Dtype* sample_mat, Dtype* dist, Dtype eps) {
    uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num && idy < column_sample_num) {
        const Dtype* bottom_row0 = bottom_data0 + idx * dim;
        const Dtype* bottom_row1 = bottom_data1 + int(sample_mat[idx * num + idy] - 1) * dim;
        Dtype norm0 = 0, norm1 = 0, dotproduct = 0;
        for (int j = 0; j < dim; ++j) {
            norm0 += bottom_row0[j] * bottom_row0[j];
            norm1 += bottom_row1[j] * bottom_row1[j];
            dotproduct += bottom_row0[j] * bottom_row1[j];
        }
        dist[idx * column_sample_num + idy] = 0.5 + 0.5 * dotproduct / (sqrt(norm0 * norm1) + eps);
    }
}

template <typename Dtype>
__global__ void L2PropagateDownKernel(const int num, const int dim, const int column_sample_num, const Dtype* bottom_data0, const Dtype* bottom_data1, const Dtype* top_diff, const Dtype* sample_mat_array0, Dtype* derivative_mat_array0, Dtype* derivative_mat_array1) {
    uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num && idy < column_sample_num) {
        const Dtype* bottom_row0 = bottom_data0 + idx * dim;
        const Dtype* bottom_row1 = bottom_data1 + int(sample_mat_array0[idx * num + idy] - 1) * dim;
        Dtype* derivative_mat_array0_row = derivative_mat_array0 + (idx * column_sample_num + idy) * dim;
        Dtype* derivative_mat_array1_row = derivative_mat_array1 + (idx * column_sample_num + idy) * dim;
        Dtype pmult = -2 * top_diff[idx * column_sample_num + idy];
        for (int j = 0; j < dim; ++j) {
            derivative_mat_array0_row[j] = pmult * (bottom_row0[j] - bottom_row1[j]);
            derivative_mat_array1_row[j] = pmult * (bottom_row1[j] - bottom_row0[j]);
        }
    }
}

template <typename Dtype>
__global__ void SumUpDerivativeKernel0(const int num, const int dim, const int column_sample_num, Dtype* derivative_mat_array0, Dtype* bottom_diff0) {
    uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num && idy < dim) {
        for (int j = 0; j < column_sample_num; j++) {
            bottom_diff0[idx * dim + idy] += derivative_mat_array0[idx * dim * column_sample_num + j * dim + idy];
        }
    }
}

template <typename Dtype>
__global__ void SumUpDerivativeKernel1(const int num, const int dim, const int column_sample_num, const Dtype* sample_mat_array1, Dtype* derivative_mat_array1, Dtype* bottom_diff1) {
    uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num && idy < dim) {
        for (int j = 0; (j < num) && (sample_mat_array1[idx * num + j] > 0); j++) {
            bottom_diff1[idx * dim + idy] += derivative_mat_array1[int(sample_mat_array1[idx * num + j] - 1) * dim + idy];
        }
    }
}

template <typename Dtype>
__global__ void CosinePropagateDownKernel(const int num, const int dim, const int column_sample_num, const Dtype* bottom_data0, const Dtype* bottom_data1, const Dtype* top_diff, const Dtype* sample_mat_array0, Dtype* derivative_mat_array0, Dtype* derivative_mat_array1, Dtype eps) {
    uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;
    uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < num && idy < column_sample_num) {
        const Dtype* bottom_row0 = bottom_data0 + idx * dim;
        const Dtype* bottom_row1 = bottom_data1 + int(sample_mat_array0[idx * num + idy] - 1) * dim;
        Dtype* derivative_mat_array0_row = derivative_mat_array0 + (idx * column_sample_num + idy) * dim;
        Dtype* derivative_mat_array1_row = derivative_mat_array1 + (idx * column_sample_num + idy) * dim;
        Dtype norm0 = 0, norm1 = 0, dotproduct = 0;
        for (int j = 0; j < dim; ++j) {
            norm0 += bottom_row0[j] * bottom_row0[j];
            norm1 += bottom_row1[j] * bottom_row1[j];
            dotproduct += bottom_row0[j] * bottom_row1[j];
        }
        norm0 = sqrt(norm0);
        norm1 = sqrt(norm1);
        Dtype tmp = norm0 * norm1;
        Dtype alpha = 0.5 / (tmp > eps ? tmp : eps);
        Dtype tmp0 = pow(norm0, Dtype(3.0)) * norm1;
        Dtype tmp1 = pow(norm1, Dtype(3.0)) * norm0;
        Dtype beta0 = -0.5 * dotproduct / (tmp0 > eps ? tmp0 : eps);
        Dtype beta1 = -0.5 * dotproduct / (tmp1 > eps ? tmp1 : eps);
        Dtype grad = top_diff[idx * column_sample_num + idy];
        for (int j = 0; j < dim; ++j) {
            derivative_mat_array0_row[j] = grad * (alpha * bottom_row1[j] + beta0 * bottom_row0[j]);
            derivative_mat_array1_row[j] = grad * (alpha * bottom_row0[j] + beta1 * bottom_row1[j]);
        }
    }
}

template <typename Dtype>
void L2JointSimilarityLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
    const Dtype* bottom_data0 = bottom[0]->gpu_data();
    const Dtype* bottom_data1 = pairwise_ ? bottom[1]->gpu_data() : bottom[0]->gpu_data();
    Dtype* dist = top[0]->mutable_gpu_data();
    int num = bottom[0]->num();
    int dim = bottom[0]->count() / bottom[0]->num();
    Sample_pairs(num, dim);
    Dtype* sample_mat_array0_gpu = sample_mat0_.mutable_gpu_data();
    Dtype* sample_mat_array1_gpu = sample_mat1_.mutable_gpu_data();
    
    dim3 thread_tail(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
    dim3 block_tail((num + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (column_sample_num_ + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
    L2DistKernel<Dtype> << <block_tail, thread_tail>> >(num, dim, column_sample_num_, bottom_data0, bottom_data1, sample_mat_array0_gpu, dist);

}

template <typename Dtype>
void CosineJointSimilarityLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
    const Dtype* bottom_data0 = bottom[0]->gpu_data();
    const Dtype* bottom_data1 = pairwise_ ? bottom[1]->gpu_data() : bottom[0]->gpu_data();
    Dtype* dist = top[0]->mutable_gpu_data();
    int num = bottom[0]->num();
    int dim = bottom[0]->count() / bottom[0]->num();
    Sample_pairs(num, dim);
    Dtype* sample_mat_array0_gpu = sample_mat0_.mutable_gpu_data();
    Dtype* sample_mat_array1_gpu = sample_mat1_.mutable_gpu_data();
   
    dim3 thread_tail(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
    dim3 block_tail((num + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (column_sample_num_ + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
    CosineDistKernel<Dtype> << <block_tail, thread_tail >> >(num, dim, column_sample_num_, bottom_data0, bottom_data1, sample_mat_array0_gpu, dist, eps_);
}

template <typename Dtype>
void L2JointSimilarityLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    const Dtype* bottom_data0 = bottom[0]->gpu_data();
    const Dtype* bottom_data1 = pairwise_ ? bottom[1]->gpu_data() : bottom[0]->gpu_data();
    const Dtype* sample_mat_array0 = sample_mat0_.gpu_data();
    const Dtype* sample_mat_array1 = sample_mat1_.gpu_data();
    const Dtype* top_diff = top[0]->gpu_diff();
    Dtype* derivative_mat_array0 = derivative_mat0_.mutable_gpu_data();
    Dtype* derivative_mat_array1 = derivative_mat1_.mutable_gpu_data();
    Dtype* bottom_diff0 = propagate_down[0] ? bottom[0]->mutable_gpu_diff() : nullptr;
    Dtype* bottom_diff1 = propagate_down[1] ? (pairwise_ ? bottom[1]->mutable_gpu_diff() : bottom[0]->mutable_gpu_diff()) : nullptr;
    if (propagate_down[0])
        caffe_gpu_set(bottom[0]->count(), Dtype(0), bottom_diff0);
    if (propagate_down[1] && pairwise_)
        caffe_gpu_set(bottom[1]->count(), Dtype(0), bottom_diff1);
    int num = bottom[0]->num();
    int dim = bottom[0]->count() / bottom[0]->num();
    // Only compute matched pair
    dim3 thread_tail(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
    dim3 block_tail((num + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (column_sample_num_ + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
    L2PropagateDownKernel<< <block_tail, thread_tail >> >(num, dim, column_sample_num_, bottom_data0, bottom_data1, top_diff, sample_mat_array0, derivative_mat_array0, derivative_mat_array1);
    dim3 thread_tail_sumup(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
    dim3 block_tail_sumup((num + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (dim + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
    if (propagate_down[0]) {
        SumUpDerivativeKernel0<< <block_tail_sumup, thread_tail_sumup >> >(num, dim, column_sample_num_, derivative_mat_array0, bottom_diff0);
    }
    if (propagate_down[1]) {
        SumUpDerivativeKernel1 << <block_tail_sumup, thread_tail_sumup >> >(num, dim, column_sample_num_, sample_mat_array1, derivative_mat_array1, bottom_diff1);
    }
}

template <typename Dtype>
void CosineJointSimilarityLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    const Dtype* bottom_data0 = bottom[0]->gpu_data();
    const Dtype* bottom_data1 = pairwise_ ? bottom[1]->gpu_data() : bottom[0]->gpu_data();
    const Dtype* sample_mat_array0 = sample_mat0_.gpu_data();
    const Dtype* sample_mat_array1 = sample_mat1_.gpu_data();
    const Dtype* top_diff = top[0]->gpu_diff();
    Dtype* derivative_mat_array0 = derivative_mat0_.mutable_gpu_data();
    Dtype* derivative_mat_array1 = derivative_mat1_.mutable_gpu_data();
    Dtype* bottom_diff0 = propagate_down[0] ? bottom[0]->mutable_gpu_diff() : nullptr;
    Dtype* bottom_diff1 = propagate_down[1] ? (pairwise_ ? bottom[1]->mutable_gpu_diff() : bottom[0]->mutable_gpu_diff()) : nullptr;
    if (propagate_down[0])
        caffe_gpu_set(bottom[0]->count(), Dtype(0), bottom_diff0);
    if (propagate_down[1] && pairwise_)
        caffe_gpu_set(bottom[1]->count(), Dtype(0), bottom_diff1);
    int num = bottom[0]->num();
    int dim = bottom[0]->count() / bottom[0]->num();
    // Only compute matched pair
    dim3 thread_tail(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
    dim3 block_tail((num + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (column_sample_num_ + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
    CosinePropagateDownKernel << <block_tail, thread_tail >> >(num, dim, column_sample_num_, bottom_data0, bottom_data1, top_diff, sample_mat_array0, derivative_mat_array0, derivative_mat_array1, eps_);
    dim3 thread_tail_sumup(CAFFE_THREAD_PER_BLOCK_DIM, CAFFE_THREAD_PER_BLOCK_DIM);
    dim3 block_tail_sumup((num + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM, (dim + CAFFE_THREAD_PER_BLOCK_DIM - 1) / CAFFE_THREAD_PER_BLOCK_DIM);
    if (propagate_down[0]) {
        SumUpDerivativeKernel0 << <block_tail_sumup, thread_tail_sumup >> >(num, dim, column_sample_num_, derivative_mat_array0, bottom_diff0);
    }
    if (propagate_down[1]) {
        SumUpDerivativeKernel1 << <block_tail_sumup, thread_tail_sumup >> >(num, dim, column_sample_num_, sample_mat_array1, derivative_mat_array1, bottom_diff1);
    }
}

INSTANTIATE_LAYER_GPU_FUNCS(L2JointSimilarityLayer);
INSTANTIATE_LAYER_GPU_FUNCS(CosineJointSimilarityLayer);

}  // namespace caffe
