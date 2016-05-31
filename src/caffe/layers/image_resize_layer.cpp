#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>

#include "caffe/layers/image_resize_layer.hpp"
#include <ppl.h>
using namespace concurrency;

namespace caffe{

    template<typename Dtype>
    void ImageResizeLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        LOG(INFO) << "Setting up ImageResizeLayer";
        ImageResizeParameter image_resize_param = this->layer_param_.image_resize_param();
        new_height_ = image_resize_param.new_height();
        new_width_ = image_resize_param.new_width();
    }

    template<typename Dtype>
    void ImageResizeLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        CHECK(top.size() == 1)
            << "There should be only one output for ImageResizeLayer.";
        CHECK(bottom.size() == 1)
            << "There should be only one output for ImageResizeLayer.";
        CHECK(bottom[0]->shape().size() == 4)
            << "Image Blob should have 4 dims.";
        CHECK(bottom[0]->channels() == 3)
            << "Image should have three channels.";
        top[0]->Reshape(bottom[0]->num(), bottom[0]->channels(),
            new_height_, new_width_);
        ori_image_holder_.Reshape(1, bottom[0]->channels(),
            bottom[0]->height(), bottom[0]->width());
        new_image_holder_.Reshape(1, bottom[0]->channels(),
            new_height_, new_width_);
    }

    template<typename Dtype>
    void ImageResizeLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        Dtype* bottom_data = bottom[0]->mutable_cpu_data();
        Dtype* top_data = top[0]->mutable_cpu_data();
        Dtype* ori_image_data = ori_image_holder_.mutable_cpu_data();
        Dtype* new_image_data = new_image_holder_.mutable_cpu_data();
        const int batch_size = bottom[0]->num();
        const int channel = bottom[0]->channels();
        const int old_height = bottom[0]->height();
        const int old_width = bottom[0]->width();
        const int old_image_size = bottom[0]->count(1);
        const int new_image_size = top[0]->count(1);
        int input_offset = 0;
        int output_offset = 0;
        for (int i = 0; i < batch_size; i++, input_offset += old_image_size, output_offset += new_image_size) {
            for (int h = 0; h < old_height; h++) {
                for (int w = 0; w < old_width; w++) {
                    for (int c = 0; c < channel; c++) {
                        ori_image_data[c + w * channel + h * channel * old_width] = bottom_data[input_offset + c * old_height * old_width + h * old_width + w];
                    }
                }
            }
            cv::Mat old_img = cv::Mat(cv::Size(old_height, old_width), CV_32FC3, ori_image_data);
            cv::Mat new_img = cv::Mat(cv::Size(new_height_, new_width_), CV_32FC3, new_image_data);
            cv::resize(old_img, new_img, new_img.size());
            for (int h = 0; h < new_height_; h++) {
                for (int w = 0; w < new_width_; w++) {
                    for (int c = 0; c < channel; c++) {
                        top_data[output_offset + c * new_height_ * new_width_ + h * new_width_ + w] = new_image_data[c + w * channel + h * channel * new_width_];
                    }
                }
            }
        }
    }

    template<typename Dtype>
    void ImageResizeLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        Forward_cpu(bottom, top);
    }

    INSTANTIATE_CLASS(ImageResizeLayer);
    REGISTER_LAYER_CLASS(ImageResize);
}