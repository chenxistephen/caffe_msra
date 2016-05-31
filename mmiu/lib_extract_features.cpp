#define __CAFFE_DLL_EXPORT
#include "lib_extract_features.h"
#include <stdio.h>  
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>  

#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <stdint.h>
#include "boost/algorithm/string.hpp"
#include "google/protobuf/text_format.h"
#include "leveldb/db.h"
#include "leveldb/write_batch.h"

#include "caffe/blob.hpp"
#include "caffe/common.hpp"
#include "caffe/net.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/util/io.hpp"
#include "caffe/layers/image_data_layer.hpp"
#include "caffe/layers/unified_data_layer.hpp"



using namespace caffe;


class CaffeFeatureExtractor: public IFeatureExtractor
{
public:
    bool LoadModel(const char* modelFilename, const char* prototxt);
    void ExtractCaffeFeature(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature);
    void ExtractCaffeImageFeature(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature);
    void ExtractCaffeTextFeature(unsigned char* pTextBuffer, const int length, const char* layerName, float* outFeature);
    void ExtractCaffeVectorFeature(unsigned char* pVectorBuffer, const int length, const char* layerName, float* outFeature);
    size_t GetFeatureDimension(const char* layerName);
    void ReadBgrBufferToDatum(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, Datum* datum);
    void Destroy();

private:
    boost::shared_ptr<Net<float> > feature_extraction_net;
    boost::shared_ptr<ImageDataLayer<float> > imageDataLayer;
};



bool CaffeFeatureExtractor::LoadModel(const char* modelFilename, const char* prototxt)
{
    openblas_set_num_threads(1);
    Caffe::set_mode(Caffe::CPU);
    
    string pretrained_binary_proto(modelFilename);
    string feature_extraction_proto(prototxt);
    feature_extraction_net.reset(new Net<float>(feature_extraction_proto, caffe::TEST));
    feature_extraction_net->CopyTrainedLayersFrom(pretrained_binary_proto);

    boost::shared_ptr<Layer<float> > dataLayer = feature_extraction_net->layers()[0];
    imageDataLayer = boost::dynamic_pointer_cast<ImageDataLayer<float>>(dataLayer);

    return true;
}

void CaffeFeatureExtractor::ExtractCaffeFeature(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature)
{
    ExtractCaffeImageFeature(pImageBgrBuffer, width, height,stride, layerName, outFeature);
}


void CaffeFeatureExtractor::ExtractCaffeImageFeature(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature)
{
    Caffe::set_mode(Caffe::CPU);

    Datum datum;
    ReadBgrBufferToDatum(pImageBgrBuffer, width, height, stride, &datum);
    
    boost::shared_ptr<Layer<float> > dataLayer = feature_extraction_net->layers()[0];
    boost::shared_ptr<ImageDataLayer<float> > imageDataLayer = boost::dynamic_pointer_cast<ImageDataLayer<float>>(dataLayer);
    
    // Get a batch from the blocking queue
    Batch<float>* batch = imageDataLayer -> GetOneBatch();

    float* top_data = batch->data_.mutable_cpu_data();
    float* top_label = batch->label_.mutable_cpu_data();
    imageDataLayer->transformed_data().set_cpu_data(top_data);
    imageDataLayer->data_transformer()->Transform(datum, &imageDataLayer->transformed_data());
    top_label[0] = datum.label();
    
    vector<Blob<float>*>* top = const_cast<vector<Blob<float>*>*>(&feature_extraction_net->top_vecs()[0]);
    caffe_copy(batch->data_.count(), batch->data_.cpu_data(), (*top)[0]->mutable_cpu_data());
    caffe_copy(batch->label_.count(), batch->label_.cpu_data(), (*top)[1]->mutable_cpu_data());

    feature_extraction_net->ForwardFromTo(1, feature_extraction_net->layers().size()-1);
    const boost::shared_ptr<Blob<float> > feature_blob = feature_extraction_net->blob_by_name(layerName);

    float* feature_blob_data = feature_blob->mutable_cpu_data() + feature_blob->offset(0);
    int dim_features = feature_blob->count();
    for (int d = 0; d < dim_features; ++d) {
        outFeature[d] = feature_blob_data[d];
    }        
}

void CaffeFeatureExtractor::ExtractCaffeTextFeature(unsigned char* pTextBuffer, const int length, const char* layerName, float* outFeature) {
    Caffe::set_mode(Caffe::CPU);

    vector<string> text_data;
    text_data.push_back(string(pTextBuffer, pTextBuffer + length));
    boost::shared_ptr<Layer<float> > dataLayer = feature_extraction_net->layers()[0];
    boost::shared_ptr<UnifiedDataLayer<float> > unifiedDataLayer = boost::dynamic_pointer_cast<UnifiedDataLayer<float>>(dataLayer);

    // Get a batch from the blocking queue
    UnifiedBatch<float>* batch = unifiedDataLayer->GetOneBatch();
    SparseBlob<float>& blob = *(batch->data_[0]->SparseBlob());
    unifiedDataLayer->data_loader_vec()[0].Load(text_data, *batch->data_[0]);
    vector<Blob<float>*>* top = const_cast<vector<Blob<float>*>*>(&feature_extraction_net->top_vecs()[0]);
    caffe_copy(blob.count(), blob.cpu_data(), (*top)[0]->SparseBlob()->mutable_cpu_data());
    feature_extraction_net->ForwardFromTo(1, feature_extraction_net->layers().size() - 1);
    const boost::shared_ptr<Blob<float> > feature_blob = feature_extraction_net->blob_by_name(layerName);
    float* feature_blob_data = feature_blob->mutable_cpu_data() + feature_blob->offset(0);
    int dim_features = feature_blob->count();
    for (int d = 0; d < dim_features; ++d) {
        outFeature[d] = feature_blob_data[d];
    }
}

void CaffeFeatureExtractor::ExtractCaffeVectorFeature(unsigned char* pVectorBuffer, const int length, const char* layerName, float* outFeature) {
    Caffe::set_mode(Caffe::CPU);
    
    boost::shared_ptr<Layer<float> > dataLayer = feature_extraction_net->layers()[0];
    boost::shared_ptr<UnifiedDataLayer<float> > unifiedDataLayer = boost::dynamic_pointer_cast<UnifiedDataLayer<float>>(dataLayer);

    // Get a batch from the blocking queue
    UnifiedBatch<float>* batch = unifiedDataLayer->GetOneBatch();
    
    int dim = batch->data_[0]->shape()[1];
    CHECK(dim == length / sizeof(float)) << "feature dimension mismatch, right: " << dim << "; cur: " << length / sizeof(float);
    vector<Blob<float>*>* top = const_cast<vector<Blob<float>*>*>(&feature_extraction_net->top_vecs()[0]);
    std::memcpy((*top)[0]->mutable_cpu_data(), pVectorBuffer, length);

    feature_extraction_net->ForwardFromTo(1, feature_extraction_net->layers().size() - 1);

    const boost::shared_ptr<Blob<float> > feature_blob = feature_extraction_net->blob_by_name(layerName);
    float* feature_blob_data = feature_blob->mutable_cpu_data() + feature_blob->offset(0);
    int dim_features = feature_blob->count();
    for (int d = 0; d < dim_features; ++d) {
        outFeature[d] = feature_blob_data[d];
    }
}

size_t CaffeFeatureExtractor::GetFeatureDimension(const char* layerName)
{
    string SlayerName(layerName);
    const boost::shared_ptr<Blob<float> > feature_blob = feature_extraction_net->blob_by_name(SlayerName);
    return feature_blob->count();
}

void CaffeFeatureExtractor::ReadBgrBufferToDatum(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, Datum* datum)
{
    cv::Mat cv_img_origin(height, width, CV_8UC3, CV_RGB(1, 1, 1));
    for (int h = 0; h < height; ++h) {
        unsigned char* tmp = pImageBgrBuffer;
        for (int w = 0; w < width; ++w) {
            for (int c = 0; c < 3; ++c) {
                cv_img_origin.at<cv::Vec3b>(h, w)[c] = *(tmp + c);
            }
            tmp += 3;
        }
        pImageBgrBuffer += stride;
    }

    cv::Mat cv_img;
    int new_width = imageDataLayer->layer_param().image_data_param().new_width();
    int new_height = imageDataLayer->layer_param().image_data_param().new_height();
    cv::resize(cv_img_origin, cv_img, cv::Size(new_height, new_width));

    datum->set_channels(3);
    datum->set_height(cv_img.rows);
    datum->set_width(cv_img.cols);
    datum->set_label(0);
    datum->clear_data();
    datum->clear_float_data();
    string* datum_string = datum->mutable_data();

    for (int c = 0; c < 3; ++c) {
        for (int h = 0; h < cv_img.rows; ++h) {
            for (int w = 0; w < cv_img.cols; ++w) {
                datum_string->push_back(static_cast<char>(cv_img.at<cv::Vec3b>(h, w)[c]));
            }
        }
    }
}

void CaffeFeatureExtractor::Destroy()
{
    delete this;
}


void* CaffeCreatFeatureExtractor()
{
	return reinterpret_cast<void*>(new CaffeFeatureExtractor());
}

bool CaffeLoadModel(void* pFeatureExtractor, const char* modelFilename, const char* prototxt){
	return reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->LoadModel(modelFilename, prototxt);
}

void CaffeExtractFeature(void* pFeatureExtractor, unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature){
    reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->ExtractCaffeFeature(pImageBgrBuffer, width, height, stride, layerName, outFeature);
}

void CaffeExtractImageFeature(void* pFeatureExtractor, unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature){
    reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->ExtractCaffeImageFeature(pImageBgrBuffer, width, height, stride, layerName, outFeature);
}

void CaffeExtractTextFeature(void* pFeatureExtractor, unsigned char* pTextBuffer, const int length, const char* layerName, float* outFeature){
    reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->ExtractCaffeTextFeature(pTextBuffer, length, layerName, outFeature);
}

void CaffeExtractVectorFeature(void* pFeatureExtractor, unsigned char* pVectorBuffer, const int length, const char* layerName, float* outFeature){
    reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->ExtractCaffeVectorFeature(pVectorBuffer, length, layerName, outFeature);
}

size_t CaffeGetFeatureDimension(void* pFeatureExtractor, const char* layerName){
    return reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->GetFeatureDimension(layerName);
}

void CaffeDestroy(void* pFeatureExtractor){
	reinterpret_cast<IFeatureExtractor*>(pFeatureExtractor)->Destroy();
}