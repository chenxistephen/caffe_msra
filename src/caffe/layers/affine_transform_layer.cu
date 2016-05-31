#include <algorithm>
#include <cfloat>
#include <vector>

#include "caffe/layer.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/layers/affine_transform_layer.hpp"

namespace caffe {

	__device__ bool is_point_in_region(int x, int y, int width_, int height_) {
		return (x >= 0 && x < width_ && y >= 0 && y < height_);
	}

	template <typename Dtype>
	__global__ void AffineTransformForward(const int count, const Dtype* bottom_data, const Dtype* transform_matrix,
		const int num, const int channels, const int height_,
		const int width_, const int output_height_, const int output_width_,
		Dtype* top_data) {
		CUDA_KERNEL_LOOP(index, count) {
			int w = index % output_width_;
			int h = (index / output_width_) % output_height_;
			int c = (index / output_width_ / output_height_) % channels;
			int n = index / output_width_ / output_height_ / channels;
			
			Dtype out_y = -1 + (Dtype)h / (output_height_ - 1) * 2;
			Dtype out_x = -1 + (Dtype)w / (output_width_ - 1) * 2;
			Dtype source_norm_x = out_x * transform_matrix[n * 6 + 0] + out_y * transform_matrix[n * 6 + 1] + transform_matrix[n * 6 + 2]; //normalized grid, -1 to 1
			Dtype source_norm_y = out_x * transform_matrix[n * 6 + 3] + out_y * transform_matrix[n * 6 + 4] + transform_matrix[n * 6 + 5]; //normalized grid, -1 to 1
			Dtype sorce_x = (source_norm_x + 1) * (width_ - 1) / 2;
			Dtype sorce_y = (source_norm_y + 1) * (height_ - 1) / 2;

			int yInTopLeft, xInTopLeft;
			Dtype yWeightTopLeft, xWeightTopLeft;
			xInTopLeft = int(sorce_x);
			yInTopLeft = int(sorce_y);
			xWeightTopLeft = 1 - (sorce_x - xInTopLeft);
			yWeightTopLeft = 1 - (sorce_y - yInTopLeft);

			// Check if the source point in the region
			bool topLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft, width_, height_);
			bool topRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft, width_, height_);
			bool bottomLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft + 1, width_, height_);
			bool bottomRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft + 1, width_, height_);

			Dtype inTopLeft = 0;
			Dtype inTopRight = 0;
			Dtype inBottomLeft = 0;
			Dtype inBottomRight = 0;

			if (topLeftIsIn) inTopLeft = *(bottom_data + ((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft);
			if (topRightIsIn) inTopRight = *(bottom_data + ((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft + 1); 
			if (bottomLeftIsIn) inBottomLeft = *(bottom_data + ((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft); 
			if (bottomRightIsIn) inBottomRight = *(bottom_data + ((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft + 1);

			Dtype v = xWeightTopLeft * yWeightTopLeft * inTopLeft
				+ (1 - xWeightTopLeft) * yWeightTopLeft * inTopRight
				+ xWeightTopLeft * (1 - yWeightTopLeft) * inBottomLeft
				+ (1 - xWeightTopLeft) * (1 - yWeightTopLeft) * inBottomRight;

			*(top_data + ((n * channels + c) * output_height_ + h) * output_width_ + w) = v;
		}
	}

	template <typename Dtype>
 void AffineTransformLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
		const vector<Blob<Dtype>*>& top) {
		const Dtype* transform_matrix = bottom[1]->gpu_data();
		// Transform matrix = 
		// n * [ data1 data2 data3 ]
		//     [ data4 data5 data6 ]
		const Dtype* bottom_data = bottom[0]->gpu_data();
		Dtype* top_data = top[0]->mutable_gpu_data();
		int count = top[0]->count();
		// We'll output the mask to top[1] if it's of size >1.
		// NOLINT_NEXT_LINE(whitespace/operators)
		AffineTransformForward<Dtype> << <CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS >> >(
			count, bottom_data, transform_matrix, bottom[0]->num(), channels_,
			height_, width_, output_height_, output_width_, top_data);
		CUDA_POST_KERNEL_CHECK;
	}

 template <typename Dtype>
 __global__ void AffineTransformBackward(const int count, const Dtype* top_diff, const Dtype* bottom_data, const Dtype* transform_matrix,
     const int num, const int channels, const int height_,
     const int width_, const int output_height_, const int output_width_,
     Dtype* bottom_diff, Dtype* transform_matrix_diff);

 template <>
 __global__ void AffineTransformBackward<float>(const int count, const float* top_diff, const float* bottom_data, const float* transform_matrix,
     const int num, const int channels, const int height_,
     const int width_, const int output_height_, const int output_width_,
     float* bottom_diff, float* transform_matrix_diff) {
     CUDA_KERNEL_LOOP(index, count) {
         int w = index % output_width_;
         int h = (index / output_width_) % output_height_;
         int c = (index / output_width_ / output_height_) % channels;
         int n = index / output_width_ / output_height_ / channels;
         float gradOutValue = top_diff[((n * channels + c) * output_height_ + h) * output_width_ + w];

         float out_y = -1 + (float)h / (output_height_ - 1) * 2;
         float out_x = -1 + (float)w / (output_width_ - 1) * 2;
         float source_norm_x = out_x * transform_matrix[n * 6 + 0] + out_y * transform_matrix[n * 6 + 1] + transform_matrix[n * 6 + 2]; //normalized grid, -1 to 1
         float source_norm_y = out_x * transform_matrix[n * 6 + 3] + out_y * transform_matrix[n * 6 + 4] + transform_matrix[n * 6 + 5]; //normalized grid, -1 to 1
         float sorce_x = (source_norm_x + 1) * (width_ - 1) / 2;
         float sorce_y = (source_norm_y + 1) * (height_ - 1) / 2;

         int yInTopLeft, xInTopLeft;
         float yWeightTopLeft, xWeightTopLeft;
         xInTopLeft = int(sorce_x);
         yInTopLeft = int(sorce_y);
         xWeightTopLeft = 1 - (sorce_x - xInTopLeft);
         yWeightTopLeft = 1 - (sorce_y - yInTopLeft);

         // Check if the source point in the region
         bool topLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft, width_, height_);
         bool topRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft, width_, height_);
         bool bottomLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft + 1, width_, height_);
         bool bottomRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft + 1, width_, height_);

         float dx = 0;
         float dy = 0;
         if (topLeftIsIn)
         {
             float topLeftDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft] * gradOutValue;
             dy += -xWeightTopLeft * topLeftDotProduct;
             dx += -yWeightTopLeft * topLeftDotProduct;
             atomicAdd(bottom_diff + ((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft, xWeightTopLeft * yWeightTopLeft * gradOutValue);
         }
         if (topRightIsIn)
         {
             float topRightDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft + 1] * gradOutValue;
             dy += -(1 - xWeightTopLeft) * topRightDotProduct;
             dx += yWeightTopLeft * topRightDotProduct;
             atomicAdd(bottom_diff + ((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft + 1, (1 - xWeightTopLeft) * yWeightTopLeft * gradOutValue);
         }
         if (bottomLeftIsIn)
         {
             float bottomLeftDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft] * gradOutValue;
             dy += xWeightTopLeft * bottomLeftDotProduct;
             dx += -(1 - yWeightTopLeft) * bottomLeftDotProduct;
             atomicAdd(bottom_diff + ((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft, xWeightTopLeft * (1 - yWeightTopLeft) * gradOutValue);
         }
         if (bottomRightIsIn)
         {
             float bottomRightDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft + 1] * gradOutValue;
             dy += (1 - xWeightTopLeft) * bottomRightDotProduct;
             dx += (1 - yWeightTopLeft) * bottomRightDotProduct;
             atomicAdd(bottom_diff + ((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft + 1, (1 - xWeightTopLeft) * (1 - yWeightTopLeft) * gradOutValue);
         }
         float d_norm_y = dy * (height_ - 1) / 2;
         float d_norm_x = dx * (width_ - 1) / 2;

         atomicAdd(transform_matrix_diff + n * 6 + 0, d_norm_x * out_x);
         atomicAdd(transform_matrix_diff + n * 6 + 1, d_norm_x * d_norm_x * out_y);
         atomicAdd(transform_matrix_diff + n * 6 + 2, d_norm_x);

         atomicAdd(transform_matrix_diff + n * 6 + 3, d_norm_y * out_x);
         atomicAdd(transform_matrix_diff + n * 6 + 4, d_norm_y * out_y);
         atomicAdd(transform_matrix_diff + n * 6 + 5, d_norm_y);
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

 template <>
 __global__ void AffineTransformBackward<double>(const int count, const double* top_diff, const double* bottom_data, const double* transform_matrix,
     const int num, const int channels, const int height_,
     const int width_, const int output_height_, const int output_width_,
     double* bottom_diff, double* transform_matrix_diff) {
     CUDA_KERNEL_LOOP(index, count) {
         int w = index % output_width_;
         int h = (index / output_width_) % output_height_;
         int c = (index / output_width_ / output_height_) % channels;
         int n = index / output_width_ / output_height_ / channels;
         double gradOutValue = top_diff[((n * channels + c) * output_height_ + h) * output_width_ + w];

         double out_y = -1 + (double)h / (output_height_ - 1) * 2;
         double out_x = -1 + (double)w / (output_width_ - 1) * 2;
         double source_norm_x = out_x * transform_matrix[n * 6 + 0] + out_y * transform_matrix[n * 6 + 1] + transform_matrix[n * 6 + 2]; //normalized grid, -1 to 1
         double source_norm_y = out_x * transform_matrix[n * 6 + 3] + out_y * transform_matrix[n * 6 + 4] + transform_matrix[n * 6 + 5]; //normalized grid, -1 to 1
         double sorce_x = (source_norm_x + 1) * (width_ - 1) / 2;
         double sorce_y = (source_norm_y + 1) * (height_ - 1) / 2;

         int yInTopLeft, xInTopLeft;
         double yWeightTopLeft, xWeightTopLeft;
         xInTopLeft = int(sorce_x);
         yInTopLeft = int(sorce_y);
         xWeightTopLeft = 1 - (sorce_x - xInTopLeft);
         yWeightTopLeft = 1 - (sorce_y - yInTopLeft);

         // Check if the source point in the region
         bool topLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft, width_, height_);
         bool topRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft, width_, height_);
         bool bottomLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft + 1, width_, height_);
         bool bottomRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft + 1, width_, height_);

         double dx = 0;
         double dy = 0;
         if (topLeftIsIn)
         {
             double topLeftDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft] * gradOutValue;
             dy += -xWeightTopLeft * topLeftDotProduct;
             dx += -yWeightTopLeft * topLeftDotProduct;
             atomicAdd_double(bottom_diff + ((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft, xWeightTopLeft * yWeightTopLeft * gradOutValue);
         }
         if (topRightIsIn)
         {
             double topRightDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft + 1] * gradOutValue;
             dy += -(1 - xWeightTopLeft) * topRightDotProduct;
             dx += yWeightTopLeft * topRightDotProduct;
             atomicAdd_double(bottom_diff + ((n * channels + c) * height_ + yInTopLeft) * width_ + xInTopLeft + 1, (1 - xWeightTopLeft) * yWeightTopLeft * gradOutValue);
         }
         if (bottomLeftIsIn)
         {
             double bottomLeftDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft] * gradOutValue;
             dy += xWeightTopLeft * bottomLeftDotProduct;
             dx += -(1 - yWeightTopLeft) * bottomLeftDotProduct;
             atomicAdd_double(bottom_diff + ((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft, xWeightTopLeft * (1 - yWeightTopLeft) * gradOutValue);
         }
         if (bottomRightIsIn)
         {
             double bottomRightDotProduct = bottom_data[((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft + 1] * gradOutValue;
             dy += (1 - xWeightTopLeft) * bottomRightDotProduct;
             dx += (1 - yWeightTopLeft) * bottomRightDotProduct;
             atomicAdd_double(bottom_diff + ((n * channels + c) * height_ + yInTopLeft + 1) * width_ + xInTopLeft + 1, (1 - xWeightTopLeft) * (1 - yWeightTopLeft) * gradOutValue);
         }
         double d_norm_y = dy * (height_ - 1) / 2;
         double d_norm_x = dx * (width_ - 1) / 2;

         atomicAdd_double(transform_matrix_diff + n * 6 + 0, d_norm_x * out_x);
         atomicAdd_double(transform_matrix_diff + n * 6 + 1, d_norm_x * d_norm_x * out_y);
         atomicAdd_double(transform_matrix_diff + n * 6 + 2, d_norm_x);

         atomicAdd_double(transform_matrix_diff + n * 6 + 3, d_norm_y * out_x);
         atomicAdd_double(transform_matrix_diff + n * 6 + 4, d_norm_y * out_y);
         atomicAdd_double(transform_matrix_diff + n * 6 + 5, d_norm_y);
     }
 }


	template <typename Dtype>
	void AffineTransformLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
		const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
		if (!propagate_down[0]) {
			return;
		}

		const Dtype* input_image_data = bottom[0]->gpu_data();
		const Dtype* transform_matrix = bottom[1]->gpu_data();
		// Transform matrix = 
		// n * [ data1 data2 data3 ]
		//     [ data4 data5 data6 ]

		const Dtype* top_diff = top[0]->gpu_diff();
		Dtype* input_image_diff = bottom[0]->mutable_gpu_diff();
		Dtype* transform_matrix_diff = bottom[1]->mutable_gpu_diff();

		const int count = bottom[0]->count();
		caffe_gpu_set(bottom[0]->count(), Dtype(0.), input_image_diff);
		caffe_gpu_set(bottom[1]->count(), Dtype(0.), transform_matrix_diff);

  AffineTransformBackward<Dtype> << <CAFFE_GET_BLOCKS(count), CAFFE_CUDA_NUM_THREADS >> >(
      count, top_diff, input_image_data, transform_matrix, bottom[0]->num(), channels_,
      height_, width_, output_height_, output_width_, input_image_diff, transform_matrix_diff);
		CUDA_POST_KERNEL_CHECK;
	}


	INSTANTIATE_LAYER_GPU_FUNCS(AffineTransformLayer);


}  // namespace caffe