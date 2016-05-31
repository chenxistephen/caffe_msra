#include <algorithm>
#include <cfloat>
#include <vector>

#include "caffe/layer.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/layers/nca_loss_layer.hpp"

namespace caffe {

template <typename Dtype>
__global__ void NcaLossForward(const int num, const int num2, const int num22, int dim, const Dtype* bottom_data, Dtype* dist) {
  CUDA_KERNEL_LOOP(t, (num2 + 1) * num22) {
    int b0 = (t / num22), b1 = (t % num22);
    if (b0 > b1) {
      if (2 * b1 == num2 - 1) return;
      b0 = num2 - b0;
      b1 = num2 - 1 - b1;
    }
    b0 *= 2; b1 *= 2;
    int e0 = (b0 + 2 < num) ? (b0 + 2) : num, e1 = (b1 + 2 < num) ? (b1 + 2) : num;
    for (int k0 = b0; k0 < e0; ++k0) {
      const Dtype* bottom_row0 = bottom_data + k0 * dim;
      for (int k1 = (b1 > k0 + 1) ? b0 : (k0 + 1); k1 < e1; ++k1) {
        const Dtype* bottom_row1 = bottom_data + k1 * dim;
        Dtype sum = 0;
        for (int j = 0; j < dim; ++j) {
           Dtype d = bottom_row0[j] - bottom_row1[j];
           sum += d*d;
        }
        dist[k1 * num + k0] = dist[k0 * num + k1] = exp(-sum);
      }
    }
  }
}

template <typename Dtype>
__global__ void NcaLossForward1(const int num, const int num2, int dim, const Dtype* bottom_data, const Dtype* bottom_data1, Dtype* dist) {
  CUDA_KERNEL_LOOP(t, num2 * num2) {
    int b0 = (t / num2), b1 = (t % num2);
    b0 *= 2; b1 *= 2;
    int e0 = (b0 + 2 < num) ? (b0 + 2) : num, e1 = (b1 + 2 < num) ? (b1 + 2) : num;
    for (int k0 = b0; k0 < e0; ++k0) {
      const Dtype* bottom_row0 = bottom_data + k0 * dim;
      for (int k1 = b0; k1 < e1; ++k1) {
        const Dtype* bottom_row1 = bottom_data1 + k1 * dim;
        Dtype sum = 0;
        for (int j = 0; j < dim; ++j) {
           Dtype d = bottom_row0[j] - bottom_row1[j];
           sum += d*d;
        }
        dist[k0 * num + k1] = exp(-sum);
      }
    }
  }
}

template <typename Dtype>
void NcaLossLayer<Dtype>::CalculateDist_gpu(const Dtype* bottom_data_gpu, const Dtype* bottom_data1_gpu, int num, int dim, std::vector<Dtype>& dist, Dtype*& dist_gpu)
{
  dist.assign(num * num, Dtype(0));
  CUDA_CHECK(cudaMalloc(&dist_gpu, dist.size() * sizeof(Dtype)));
  CUDA_CHECK(cudaMemset(dist_gpu, 0, dist.size() * sizeof(Dtype)));
  if (pairwise_) {
    int num2 = (num + 1) / 2;
    NcaLossForward1<Dtype><<<CAFFE_GET_BLOCKS(num2 * num2), CAFFE_CUDA_NUM_THREADS>>>(
        num, num2, dim, bottom_data_gpu, bottom_data1_gpu, dist_gpu);
    CUDA_POST_KERNEL_CHECK;
  }
  else {
    int num2 = (num + 1) / 2;
    int num22 = (num2 + 1) / 2;
    NcaLossForward<Dtype><<<CAFFE_GET_BLOCKS((num2 + 1) * num22), CAFFE_CUDA_NUM_THREADS>>>(
        num, num2, num22, dim, bottom_data_gpu, dist_gpu);
    CUDA_POST_KERNEL_CHECK;
  }
  CUDA_CHECK(cudaMemcpy(dist.data(), dist_gpu, dist.size() * sizeof(Dtype), cudaMemcpyDeviceToHost));
}

template <typename Dtype>
void NcaLossLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  const Dtype* bottom_data_gpu = bottom[0]->gpu_data();
  const Dtype* bottom_data1_gpu = pairwise_ ? bottom[1]->gpu_data() : nullptr;
  const Dtype* bottom_label = !pairwise_ ? bottom[1]->cpu_data() : nullptr;
  int num = bottom[0]->num();
  int dim = bottom[0]->count() / bottom[0]->num();
  std::vector<Dtype> dist;
  Dtype* dist_gpu = nullptr;
  CalculateDist_gpu(bottom_data_gpu, bottom_data1_gpu, num, dim, dist, dist_gpu);
  CUDA_CHECK(cudaFree(dist_gpu));
  top[0]->mutable_cpu_data()[0] = ForwardInternal(dist.data(), bottom_label, num);
}

template <typename Dtype>
__global__ void NcaLossBackward(const int num, const int dim, const int dim8, const Dtype* dist, const Dtype* dist_grad, const Dtype* bottom_data, Dtype* bottom_diff) {
  CUDA_KERNEL_LOOP(t, num * dim8) {
    int k0 = (t / dim8), b = (t % dim8) * 8;
    int e = (b + 8 < dim) ? (b + 8) : dim;
    const Dtype* bottom_row0 = bottom_data + k0 * dim;
    Dtype* bottom_grad_row0 = bottom_diff + k0 * dim;
    for (int k1 = 0; k1 < num; ++k1) {
      if (k1 != k0) {
        const Dtype* bottom_row1 = bottom_data + k1 * dim;
        Dtype pmult = -(dist_grad[k1 * num + k0] + dist_grad[k0 * num + k1]) * dist[k0 * num + k1] * 2;
        for (int i = b; i < e; ++i) {
          bottom_grad_row0[i] += pmult * (bottom_row0[i] - bottom_row1[i]);
        }
      }
    }
  }
}

template <typename Dtype>
__global__ void NcaLossBackward1(const int num, const int dim, const int dim8, const Dtype* dist, const Dtype* dist_grad, const Dtype* bottom_data, const Dtype* bottom_data1, Dtype* bottom_diff) {
  CUDA_KERNEL_LOOP(t, num * dim8) {
    int k0 = (t / dim8), b = (t % dim8) * 8;
    int e = (b + 8 < dim) ? (b + 8) : dim;
    const Dtype* bottom_row0 = bottom_data + k0 * dim;
    Dtype* bottom_grad_row0 = bottom_diff + k0 * dim;
    for (int k1 = 0; k1 < num; ++k1) {
      const Dtype* bottom_row1 = bottom_data1 + k1 * dim;
      Dtype pmult = -dist_grad[k0 * num + k1] * dist[k0 * num + k1] * 2;
      for (int i = b; i < e; ++i) {
        bottom_grad_row0[i] += pmult * (bottom_row0[i] - bottom_row1[i]);
      }
    }
  }
}

template <typename Dtype>
__global__ void NcaLossBackward2(const int num, const int dim, const int dim8, const Dtype* dist, const Dtype* dist_grad, const Dtype* bottom_data, const Dtype* bottom_data1, Dtype* bottom_diff1) {
  CUDA_KERNEL_LOOP(t, num * dim8) {
    int k1 = (t / dim8), b = (t % dim8) * 8;
    int e = (b + 8 < dim) ? (b + 8) : dim;
    const Dtype* bottom_row1 = bottom_data1 + k1 * dim;
    Dtype* bottom_grad_row1 = bottom_diff1 + k1 * dim;
    for (int k0 = 0; k0 < num; ++k0) {
      const Dtype* bottom_row0 = bottom_data + k0 * dim;
      Dtype pmult = -dist_grad[k0 * num + k1] * dist[k0 * num + k1] * 2;
      for (int i = b; i < e; ++i) {
        bottom_grad_row1[i] += pmult * (bottom_row1[i] - bottom_row0[i]);
      }
    }
  }
}

template <typename Dtype>
void NcaLossLayer<Dtype>::Backward_gpu(
    const vector<Blob<Dtype>*>& top, const vector<bool>& propagate_down,
    const vector<Blob<Dtype>*>& bottom) {
  if (propagate_down[1] && !pairwise_) {
    LOG(FATAL) << this->type()
               << " Layer cannot backpropagate to label inputs.";
  }
  if (propagate_down[0] || propagate_down[1]) {
    const Dtype* bottom_data_gpu = bottom[0]->gpu_data();
    const Dtype* bottom_data1_gpu = pairwise_ ? bottom[1]->gpu_data() : nullptr;
    const Dtype* bottom_label = !pairwise_ ? bottom[1]->cpu_data() : nullptr;
    Dtype* bottom_diff_gpu = propagate_down[0] ? bottom[0]->mutable_gpu_diff() : nullptr;
    Dtype* bottom_diff1_gpu = propagate_down[1] ? bottom[1]->mutable_gpu_diff() : nullptr;
    int num = bottom[0]->num();
    int dim = bottom[0]->count() / bottom[0]->num();
    std::vector<Dtype> dist;
    Dtype* dist_gpu = nullptr;
    CalculateDist_gpu(bottom_data_gpu, bottom_data1_gpu, num, dim, dist, dist_gpu);
    std::vector<Dtype> dist_grad;
    BackwardInternal(dist.data(), bottom_label, num, top[0]->cpu_diff()[0], dist_grad);
    Dtype* dist_grad_gpu = nullptr;
    CUDA_CHECK(cudaMalloc(&dist_grad_gpu, dist_grad.size() * sizeof(Dtype)));
    CUDA_CHECK(cudaMemcpy(dist_grad_gpu, dist_grad.data(), dist_grad.size() * sizeof(Dtype), cudaMemcpyHostToDevice));
    if (pairwise_) {
      int dim8 = (dim + 7) / 8;
      if (propagate_down[0]) {
        CUDA_CHECK(cudaMemset(bottom_diff_gpu, 0, bottom[0]->count() * sizeof(Dtype)));
        NcaLossBackward1<Dtype><<<CAFFE_GET_BLOCKS(num * dim8), CAFFE_CUDA_NUM_THREADS>>>(
            num, dim, dim8, dist_gpu, dist_grad_gpu, bottom_data_gpu, bottom_data1_gpu, bottom_diff_gpu);
        CUDA_POST_KERNEL_CHECK;
      }
      if (propagate_down[1]) {
        CUDA_CHECK(cudaMemset(bottom_diff1_gpu, 0, bottom[1]->count() * sizeof(Dtype)));
        NcaLossBackward2<Dtype><<<CAFFE_GET_BLOCKS(num * dim8), CAFFE_CUDA_NUM_THREADS>>>(
            num, dim, dim8, dist_gpu, dist_grad_gpu, bottom_data_gpu, bottom_data1_gpu, bottom_diff1_gpu);
        CUDA_POST_KERNEL_CHECK;
      }
    }
    else {
      int dim8 = (dim + 7) / 8;
      CUDA_CHECK(cudaMemset(bottom_diff_gpu, 0, bottom[0]->count() * sizeof(Dtype)));
      NcaLossBackward<Dtype><<<CAFFE_GET_BLOCKS(num * dim8), CAFFE_CUDA_NUM_THREADS>>>(
          num, dim, dim8, dist_gpu, dist_grad_gpu, bottom_data_gpu, bottom_diff_gpu);
      CUDA_POST_KERNEL_CHECK;
    }
    CUDA_CHECK(cudaFree(dist_gpu));
    CUDA_CHECK(cudaFree(dist_grad_gpu));
  }
}

INSTANTIATE_LAYER_GPU_FUNCS(NcaLossLayer);

}  // namespace caffe
