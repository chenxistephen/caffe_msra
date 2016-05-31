#include "caffe/image_data_transformer.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>
#include <opencv2/imgproc/imgproc.hpp>
#include "IUNativeExecutor.h"
#include "caffe/util/rng.hpp"
#include "caffe/util/io.hpp"
#include "ppl.h"
using namespace concurrency;
#include <fstream>

namespace caffe {

    void write_binary_file(const char* file_name, const uint8_t* data, int len) {
        std::ofstream out(file_name, ios::binary);
        out.write((char*)data, len);
        out.close();
    }

    std::string& replace_string(std::string &s,
        const std::string & olds,
        const std::string & news)
    {
        auto it = s.find(olds);
        while (it != std::string::npos) {
            s = s.replace(it, olds.length(), news);
            it = s.find(olds);
        }
        return s;
    }

    template<class Dtype>
    shared_ptr<IImageDataTransformer<Dtype>> IImageDataTransformer<Dtype>::CreateInstanceAndInit(const string& type, const ImageDataTransformParameter& param, Phase phase){
        shared_ptr<IImageDataTransformer<Dtype>> obj;
        if (type == "caffe")
            obj.reset(new CaffeImageDataTransformer<Dtype>());
        if (type == "iu")
            obj.reset(new IUImageDataTransformer<Dtype>());
        CHECK(obj.get() != nullptr) << "non-supported algorithm: " << type;
        obj->Init(param, phase);
        return obj;
    }

    template<class Dtype>
    void CaffeImageDataTransformer<Dtype>::Init(const ImageDataTransformParameter& param, Phase phase) {
        param_ = param;
        phase_ = phase;
        data_transformer_.reset(
            new DataTransformer<Dtype>(param.transform_param(), this->phase_));
    }

    template<class Dtype>
    void CaffeImageDataTransformer<Dtype>::Transform(const vector<ImageData> &image_data, Blob<Dtype>& out_data) const {
        //for each(auto& d in image_data) CHECK(d.bbox.IsNull()) << "caffe transformer doesn't support bbox";

        bool needs_resize = param_.transform_param().has_resize_width() && param_.transform_param().has_resize_height();
        const int crop_size = param_.transform_param().crop_size();

        CHECK(!(param_.transform_param().has_resize_width() ^ param_.transform_param().has_resize_height())) << "resize width and height need to be specified simultaneously";

        data_transformer_->InitRand();
        cv::Mat cv_image;
        Blob<Dtype> transformed_data;
        for (size_t i = 0; i < image_data.size(); ++i) {
            cv_image = cv::imdecode(image_data[i].data, CV_LOAD_IMAGE_COLOR);
            if (i == 0) {
                const int channels = cv_image.channels();
                const int height = cv_image.rows;
                const int width = cv_image.cols;
                if (crop_size > 0) {
                    out_data.Reshape((int)image_data.size(), channels, crop_size, crop_size);
                    transformed_data.Reshape(1, channels, crop_size, crop_size);
                }
                else if (needs_resize) {
                    out_data.Reshape((int)image_data.size(), channels, param_.transform_param().resize_height(), param_.transform_param().resize_width());
                    transformed_data.Reshape(1, channels, param_.transform_param().resize_height(), param_.transform_param().resize_width());
                }
                else{
                    out_data.Reshape((int)image_data.size(), channels, height, width);
                    transformed_data.Reshape(1, channels, height, width);
                }
                
            }
            CHECK(cv_image.rows * cv_image.cols > 0) << "incorrect image: " << image_data[i].key;
            int offset = out_data.offset(i);
            transformed_data.set_cpu_data(out_data.mutable_cpu_data() + offset);

            //LOG(ERROR) << "bbox: " << image_data[i].bbox.left << "," << image_data[i].bbox.top << "," << image_data[i].bbox.width << "," << image_data[i].bbox.height << "; valid: " << !image_data[i].bbox.IsNull();

            if (!image_data[i].bbox.IsNull()) {
                CHECK(image_data[i].bbox.width >= 1 && image_data[i].bbox.height >= 1 && image_data[i].bbox.left + image_data[i].bbox.width <= cv_image.cols && image_data[i].bbox.top + image_data[i].bbox.height <= cv_image.rows)
                    << "incorrect bbox params: " << image_data[i].bbox.left << "," << image_data[i].bbox.top << "," << image_data[i].bbox.width << "," << image_data[i].bbox.height;
                cv::Rect roi(image_data[i].bbox.left, image_data[i].bbox.top, image_data[i].bbox.width, image_data[i].bbox.height);
                cv::Mat cv_cropped_img = cv_image(roi);

                CHECK(needs_resize) << "needs to specify resize for inputs with bbox";
                cv::Size cv_new_size(param_.transform_param().resize_width(), param_.transform_param().resize_height());
                cv::resize(cv_cropped_img, cv_image,
                    cv_new_size, 0, 0, cv::INTER_LINEAR);

                //string image_key = image_data[i].key;   cv::imwrite(replace_string(replace_string(image_key, "\\", "_"), ":", "_") + ".cv.crop.bmp", cv_cropped_img);
            }
            else {
                if (needs_resize) {
                    cv::Size cv_new_size(param_.transform_param().resize_width(), param_.transform_param().resize_height());
                    cv::resize(cv_image, cv_image,
                        cv_new_size, 0, 0, cv::INTER_LINEAR);
                }
            }
            //string image_key = image_data[i].key;   cv::imwrite(replace_string(replace_string(image_key, "\\", "_"), ":", "_") + ".cv.bmp", cv_image);
            data_transformer_->Transform(cv_image, &transformed_data);
        }
    }

    template<class Dtype>
    void IUImageDataTransformer<Dtype>::Init(const ImageDataTransformParameter& param, Phase phase) {
        param_ = param;
        phase_ = phase;

        CHECK(!param_.transform_param().has_resize_width() && !param_.transform_param().has_resize_height()) << "IUImageDataTransformer doesn't support resize size in prototxt, you should specify it in the IU pipeline file";

        // check if we want to use mean_file
        if (param_.transform_param().has_mean_file()) {
            CHECK_EQ(param_.transform_param().mean_value_size(), 0) <<
                "Cannot specify mean_file and mean_value at the same time";
            const string& mean_file = param_.transform_param().mean_file();
            LOG(INFO) << "Loading mean file from: " << mean_file;
            BlobProto blob_proto;
            ReadProtoFromBinaryFileOrDie(mean_file.c_str(), &blob_proto);
            data_mean_.FromProto(blob_proto);
        }
        // check if we want to use mean_value
        if (param_.transform_param().mean_value_size() > 0) {
            CHECK(param_.transform_param().has_mean_file() == false) <<
                "Cannot specify mean_file and mean_value at the same time";
            for (int c = 0; c < param_.transform_param().mean_value_size(); ++c) {
                mean_values_.push_back(param_.transform_param().mean_value(c));
            }
        }

        iu_executor_.reset(new IUNativeExecutorWrapper::IUNativeExecutor);
        const char* input_names[] { "ImageStream", "Bbox"};
        const char* output_names[] { "Thumb" };
        iu_executor_->Initialize(param_.iu_pipeline().c_str(), input_names, 2, output_names, 1, ".");
        InitRand();
    }

    template <typename Dtype>
    int IUImageDataTransformer<Dtype>::Rand(int n) const {
        CHECK(rng_);
        CHECK_GT(n, 0);
        cs_.lock();
        caffe::rng_t* rng =
            static_cast<caffe::rng_t*>(rng_->generator());
        int v = ((*rng)() % n);
        cs_.unlock();
        return v;
    }

    template<class Dtype>
    void IUImageDataTransformer<Dtype>::Transform(const cv::Mat& cv_img, Blob<Dtype>& transformed_blob) const {
        const int img_channels = cv_img.channels();
        const int img_height = cv_img.rows;
        const int img_width = cv_img.cols;

        const int channels = transformed_blob.channels();
        const int height = transformed_blob.height();
        const int width = transformed_blob.width();
        const int num = transformed_blob.num();

        CHECK_EQ(channels, img_channels);
        CHECK_LE(height, img_height);
        CHECK_LE(width, img_width);
        CHECK_GE(num, 1);

        CHECK(cv_img.depth() == CV_8U) << "Image data type must be unsigned byte";

        const int crop_size = param_.transform_param().crop_size();
        const Dtype scale = param_.transform_param().scale();
        const bool do_mirror = param_.transform_param().mirror() && Rand(2);
        const bool has_mean_file = param_.transform_param().has_mean_file();
        const bool has_mean_values = mean_values_.size() > 0;

        CHECK_GT(img_channels, 0);
        CHECK_GE(img_height, crop_size);
        CHECK_GE(img_width, crop_size);

        const Dtype* mean = NULL;
        vector<Dtype> mean_values = mean_values_;
        if (has_mean_file) {
            CHECK_EQ(img_channels, data_mean_.channels());
            CHECK_EQ(img_height, data_mean_.height());
            CHECK_EQ(img_width, data_mean_.width());
            mean = data_mean_.cpu_data();
        }
        if (has_mean_values) {
            CHECK(mean_values_.size() == 1 || mean_values_.size() == img_channels) <<
                "Specify either 1 mean_value or as many as channels: " << img_channels;
            if (img_channels > 1 && mean_values_.size() == 1) {
                // Replicate the mean_value for simplicity
                for (int c = 1; c < img_channels; ++c) {
                    mean_values.push_back(mean_values_[0]);
                }
            }
        }

        int h_off = 0;
        int w_off = 0;
        cv::Mat cv_cropped_img = cv_img;
        if (crop_size) {
            CHECK_EQ(crop_size, height);
            CHECK_EQ(crop_size, width);
            // We only do random crop when we do training.
            if (phase_ == TRAIN) {
                h_off = Rand(img_height - crop_size + 1);
                w_off = Rand(img_width - crop_size + 1);
            }
            else {
                h_off = (img_height - crop_size) / 2;
                w_off = (img_width - crop_size) / 2;
            }
            cv::Rect roi(w_off, h_off, crop_size, crop_size);
            cv_cropped_img = cv_img(roi);
        }
        else {
            CHECK_EQ(img_height, height);
            CHECK_EQ(img_width, width);
        }

        CHECK(cv_cropped_img.data);

        Dtype* transformed_data = transformed_blob.mutable_cpu_data();
        int top_index;
        for (int h = 0; h < height; ++h) {
            const uchar* ptr = cv_cropped_img.ptr<uchar>(h);
            int img_index = 0;
            for (int w = 0; w < width; ++w) {
                for (int c = 0; c < img_channels; ++c) {
                    if (do_mirror) {
                        top_index = (c * height + h) * width + (width - 1 - w);
                    }
                    else {
                        top_index = (c * height + h) * width + w;
                    }
                    // int top_index = (c * height + h) * width + w;
                    Dtype pixel = static_cast<Dtype>(ptr[img_index++]);
                    if (has_mean_file) {
                        int mean_index = (c * img_height + h_off + h) * img_width + w_off + w;
                        transformed_data[top_index] =
                            (pixel - mean[mean_index]) * scale;
                    }
                    else {
                        if (has_mean_values) {
                            transformed_data[top_index] =
                                (pixel - mean_values[c]) * scale;
                        }
                        else {
                            transformed_data[top_index] = pixel * scale;
                        }
                    }
                }
            }
        }
    }

    template <typename Dtype>
    void IUImageDataTransformer<Dtype>::InitRand() {
        const bool needs_rand = param_.transform_param().mirror() ||
            (phase_ == TRAIN && param_.transform_param().crop_size());
        if (needs_rand) {
            const unsigned int rng_seed = caffe_rng_rand();
            rng_.reset(new Caffe::RNG(rng_seed));
        }
        else {
            rng_.reset();
        }
    }

    template<class Dtype>
    void IUImageDataTransformer<Dtype>::ProcessSingle(const ImageData &image_data, Blob<Dtype>& out_data) const {
        auto input_feature_set = new IUNativeExecutorWrapper::IUReadOnlyFeatureWrapper*[2];
        input_feature_set[0] = new IUNativeExecutorWrapper::IUReadOnlyFeatureWrapper();
        input_feature_set[0]->Set(&image_data.data[0], image_data.data.size(), image_data.data.size(), 1, 1);
        input_feature_set[1] = new IUNativeExecutorWrapper::IUReadOnlyFeatureWrapper();
        uint32_t bbox[4] {image_data.bbox.left, image_data.bbox.top, image_data.bbox.width, image_data.bbox.height};
        input_feature_set[1]->Set((const BYTE*)&bbox, sizeof(bbox), sizeof(bbox)/sizeof(bbox[0]), 1, 1);

        auto output_feature_set = new IUNativeExecutorWrapper::IUReadOnlyFeatureWrapper*[1];
        output_feature_set[0] = new IUNativeExecutorWrapper::IUReadOnlyFeatureWrapper();

        iu_executor_->Process(input_feature_set, 2, output_feature_set, 1);
        CHECK(output_feature_set[0]->GetLengthInBytes() > 0) << "IU processing falied to output the resized image: "<< image_data.key;

        std::string image_key = image_data.key;
        //write_binary_file((replace_string(replace_string(image_key, "\\", "_"), ":", "_")+".bmp").c_str(), output_feature_set[0]->GetBuffer(), output_feature_set[0]->GetLengthInBytes());

        cv::Mat cv_image = cv::imdecode(cv::_InputArray(output_feature_set[0]->GetBuffer(), output_feature_set[0]->GetLengthInBytes()), CV_LOAD_IMAGE_COLOR);

        const int channels = cv_image.channels();
        const int height = cv_image.rows;
        const int width = cv_image.cols;

        CHECK(channels == 3) << "the # channel of the image should be 3";

        const int crop_size = param_.transform_param().crop_size();
        if (out_data.count() == 0) {
            if (crop_size > 0) {
                out_data.Reshape((int)1, channels, crop_size, crop_size);
            }
            else {
                out_data.Reshape((int)1, channels, height, width);
            }
        }

        Transform(cv_image, out_data);

        delete input_feature_set[0];
        delete[] input_feature_set;
        delete output_feature_set[0];
        delete[] output_feature_set;
    }

    template<class Dtype>
    void IUImageDataTransformer<Dtype>::Transform(const vector<ImageData> &image_data, Blob<Dtype>& out_data) const {
        (const_cast<IUImageDataTransformer<Dtype>*>(this))->InitRand();
        vector<int> dim1, dim;
		
        //process the first
        {
            Blob<Dtype> transformed_data;
            ProcessSingle(image_data[0], transformed_data);
            dim1 = transformed_data.shape();
            dim = dim1;
            dim[0] = (int)image_data.size();
            out_data.Reshape(dim);
            caffe_copy(transformed_data.count(), transformed_data.cpu_data(), out_data.mutable_cpu_data() + out_data.offset(0));
        }
		
#ifndef NDEBUG
        for (size_t i = 1; i < image_data.size(); ++i)
#else
        parallel_for((size_t)1, image_data.size(), [&](size_t i)
#endif
        {
            Blob<Dtype> transformed_data;
            transformed_data.Reshape(dim1);
            int offset = out_data.offset(i);
            transformed_data.set_cpu_data(out_data.mutable_cpu_data() + offset);
            ProcessSingle(image_data[i], transformed_data);
        }
#ifdef NDEBUG
        );
#endif
    }

    INSTANTIATE_CLASS(IImageDataTransformer);
    INSTANTIATE_CLASS(CaffeImageDataTransformer);
    INSTANTIATE_CLASS(IUImageDataTransformer);
}