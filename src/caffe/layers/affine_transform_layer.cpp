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
void AffineTransformLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  AffineTransformParameter affine_transform_param = this->layer_param_.affine_transform_param();
  output_height_ = bottom[0]->height();
  output_width_ = bottom[0]->width();
  if (affine_transform_param.has_output_height())
      output_height_ = affine_transform_param.output_height();
  if (affine_transform_param.has_output_width())
      output_width_ = affine_transform_param.output_width();
}

template <typename Dtype>
void AffineTransformLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  CHECK_EQ(6 * bottom[1]->num(), bottom[1]->count()); //Affine transformation need 6 parameters
  CHECK_EQ(4, bottom[0]->num_axes()) << "Input must have 4 axes, "
      << "corresponding to (num, channels, height, width)";
  channels_ = bottom[0]->channels();
  height_ = bottom[0]->height();
  width_ = bottom[0]->width();
  top[0]->Reshape(bottom[0]->num(), channels_, output_height_,
      output_width_);
}

template <typename Dtype>
bool AffineTransformLayer<Dtype>::is_point_in_region(int x, int y) {
    return (x >= 0 && x < width_ && y >= 0 && y < height_);
}

template <typename Dtype>
void AffineTransformLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
      const vector<Blob<Dtype>*>& top) {
  const Dtype* transform_matrix = bottom[1]->cpu_data();
  // Transform matrix = 
  // n * [ data1 data2 data3 ]
  //     [ data4 data5 data6 ]
  Dtype* top_data = top[0]->mutable_cpu_data();
  for (int n = 0; n < bottom[0]->num(); ++n) {
      for (int c = 0; c < channels_; ++c) {
          for (int h = 0; h < output_height_; ++h) {
              //output grid coordinate is a normalized grid from -1 to 1
              for (int w = 0; w < output_width_; ++w) {
                  // output grid coordinate
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
                  bool topLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft);
                  bool topRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft);
                  bool bottomLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft + 1);
                  bool bottomRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft + 1);

                  Dtype inTopLeft = 0;
                  Dtype inTopRight = 0;
                  Dtype inBottomLeft = 0;
                  Dtype inBottomRight = 0;

                  if (topLeftIsIn) inTopLeft = bottom[0]->data_at(n, c, yInTopLeft, xInTopLeft);
                  if (topRightIsIn) inTopRight = bottom[0]->data_at(n, c, yInTopLeft, xInTopLeft + 1);
                  if (bottomLeftIsIn) inBottomLeft = bottom[0]->data_at(n, c, yInTopLeft + 1, xInTopLeft);
                  if (bottomRightIsIn) inBottomRight = bottom[0]->data_at(n, c, yInTopLeft + 1, xInTopLeft + 1);

                  Dtype v = xWeightTopLeft * yWeightTopLeft * inTopLeft
                      + (1 - xWeightTopLeft) * yWeightTopLeft * inTopRight
                      + xWeightTopLeft * (1 - yWeightTopLeft) * inBottomLeft
                      + (1 - xWeightTopLeft) * (1 - yWeightTopLeft) * inBottomRight;

                  top_data[top[0]->offset(n, c, h, w)] = v;
              }
          }
      }
  }
}

template <typename Dtype>
void AffineTransformLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
    const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom) {
    if (!propagate_down[0]) {
        return;
    }
    const Dtype* input_image_data = bottom[0]->cpu_data();
    const Dtype* transform_matrix = bottom[1]->cpu_data();
    // Transform matrix = 
    // n * [ data1 data2 data3 ]
    //     [ data4 data5 data6 ]

    const Dtype* top_diff = top[0]->cpu_diff();
    Dtype* input_image_diff = bottom[0]->mutable_cpu_diff();
    Dtype* transform_matrix_diff = bottom[1]->mutable_cpu_diff();

    //Clear Backward Propgationd data
    caffe_memset(bottom[0]->count() * sizeof(Dtype), 0, input_image_diff);
    caffe_memset(bottom[1]->count() * sizeof(Dtype), 0, transform_matrix_diff);

    //Backward propgation to input image

    for (int n = 0; n < bottom[0]->num(); ++n) {
        for (int c = 0; c < channels_; ++c) {
            for (int h = 0; h < output_height_; ++h) {
                //output grid coordinate is a normalized grid from -1 to 1
                for (int w = 0; w < output_width_; ++w)    {
                    Dtype gradOutValue = top[0]->diff_at(n, c, h, w);
                    // output grid coordinate
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
                    bool topLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft);
                    bool topRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft);
                    bool bottomLeftIsIn = is_point_in_region(xInTopLeft, yInTopLeft + 1);
                    bool bottomRightIsIn = is_point_in_region(xInTopLeft + 1, yInTopLeft + 1);
                    Dtype dx = 0;
                    Dtype dy = 0;
                    if (topLeftIsIn)
                    {
                        Dtype topLeftDotProduct = bottom[0]->data_at(n, c, yInTopLeft, xInTopLeft) * gradOutValue;
                        dy += -xWeightTopLeft * topLeftDotProduct;
                        dx += -yWeightTopLeft * topLeftDotProduct;
                        input_image_diff[bottom[0]->offset(n, c, yInTopLeft, xInTopLeft)] += xWeightTopLeft * yWeightTopLeft * gradOutValue;
                    }
                    if (topRightIsIn)
                    {
                        Dtype topRightDotProduct = bottom[0]->data_at(n, c, yInTopLeft, xInTopLeft + 1) * gradOutValue;
                        dy += -(1 - xWeightTopLeft) * topRightDotProduct;
                        dx += yWeightTopLeft * topRightDotProduct;
                        input_image_diff[bottom[0]->offset(n, c, yInTopLeft, xInTopLeft + 1)] += (1 - xWeightTopLeft) * yWeightTopLeft * gradOutValue;
                    }
                    if (bottomLeftIsIn)
                    {
                        Dtype bottomLeftDotProduct = bottom[0]->data_at(n, c, yInTopLeft + 1, xInTopLeft) * gradOutValue;
                        dy += xWeightTopLeft * bottomLeftDotProduct;
                        dx += -(1- yWeightTopLeft) * bottomLeftDotProduct;
                        input_image_diff[bottom[0]->offset(n, c, yInTopLeft + 1, xInTopLeft)] += xWeightTopLeft * (1 - yWeightTopLeft) * gradOutValue;
                    }
                    if (bottomRightIsIn)
                    {
                        Dtype bottomRightDotProduct = bottom[0]->data_at(n, c, yInTopLeft + 1, xInTopLeft + 1) * gradOutValue;
                        dy += (1 - xWeightTopLeft) * bottomRightDotProduct;
                        dx += (1 - yWeightTopLeft) * bottomRightDotProduct;
                        input_image_diff[bottom[0]->offset(n, c, yInTopLeft + 1, xInTopLeft + 1)] += (1 - xWeightTopLeft) * (1 - yWeightTopLeft) * gradOutValue;
                    }
                    Dtype d_norm_y = dy * (height_ - 1) / 2;
                    Dtype d_norm_x = dx * (width_ - 1) / 2;
                    transform_matrix_diff[n * 6 + 0] += d_norm_x * out_x;
                    transform_matrix_diff[n * 6 + 1] += d_norm_x * out_y;
                    transform_matrix_diff[n * 6 + 2] += d_norm_x;

                    transform_matrix_diff[n * 6 + 3] += d_norm_y * out_x;
                    transform_matrix_diff[n * 6 + 4] += d_norm_y * out_y;
                    transform_matrix_diff[n * 6 + 5] += d_norm_y;
                }
            }
        }
    }
}

#ifdef CPU_ONLY
STUB_GPU(AffineTransformLayer);
#endif

INSTANTIATE_CLASS(AffineTransformLayer);
REGISTER_LAYER_CLASS(AffineTransform);

}  // namespace caffe
