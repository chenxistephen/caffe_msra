#include <string>
#include <vector>
#include <unordered_map>

#include "caffe/common.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/blob.hpp"
#include "caffe/data_transformer.hpp"
#include "caffe/image_db.hpp"
#include "ppl.h"

namespace IUNativeExecutorWrapper {
    class IUNativeExecutor;
}

namespace caffe{

    template<class Dtype>
    class IImageDataTransformer{
    public:
        virtual ~IImageDataTransformer() {}
        virtual void Init(const ImageDataTransformParameter& param, Phase phase) = 0;
        virtual void Transform(const vector<ImageData> &image_data, Blob<Dtype>& out_data) const = 0;

    public:
        static shared_ptr<IImageDataTransformer<Dtype>> CreateInstanceAndInit(const string& type, const ImageDataTransformParameter& param, Phase phase);
    };

    template<class Dtype>
    class CaffeImageDataTransformer : public IImageDataTransformer<Dtype>{
    public:
        void Init(const ImageDataTransformParameter& param, Phase phase) override;
        void Transform(const vector<ImageData> &image_data, Blob<Dtype>& out_data) const override;

    private:
        ImageDataTransformParameter param_;
        shared_ptr<DataTransformer<Dtype>> data_transformer_;
        Phase phase_;
    };

    template<class Dtype>
    class IUImageDataTransformer : public IImageDataTransformer<Dtype>{
    public:
        void Init(const ImageDataTransformParameter& param, Phase phase) override;
        void Transform(const vector<ImageData> &image_data, Blob<Dtype>& out_data) const override;

    protected:
        void InitRand();
        int Rand(int n) const;

    private:
        void ProcessSingle(const ImageData &image_data, Blob<Dtype>& out_data) const;
        void Transform(const cv::Mat& cv_img, Blob<Dtype>& transformed_blob) const;
        mutable Concurrency::critical_section cs_;

    private:
        ImageDataTransformParameter param_;
        shared_ptr<Caffe::RNG> rng_;
        Blob<Dtype> data_mean_;
        vector<Dtype> mean_values_;
        Phase phase_;
        shared_ptr<IUNativeExecutorWrapper::IUNativeExecutor> iu_executor_;
    };
}