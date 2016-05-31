#include <fstream>  // NOLINT(readability/streams)
#include <iostream>  // NOLINT(readability/streams)
#include <string>
#include <utility>
#include <vector>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>

#include <boost/timer/timer.hpp>

#include "caffe/layers/unified_data_layer.hpp"
#include "caffe/util/benchmark.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/math_functions.hpp"
#include "caffe/util/rng.hpp"
#include "caffe/internal_thread.hpp"
#include "caffe/util/blocking_queue.hpp"
#include "caffe/common.hpp"

namespace caffe{

template <typename Dtype>
UnifiedDataLayer<Dtype>::UnifiedDataLayer(const LayerParameter& param)
  : Layer<Dtype>(param), reader_(param) {
    for (int i = 0; i < PREFETCH_COUNT; ++i) {
      prefetch_free_.push(&prefetch_[i]);
    }
  }

template <typename Dtype>
void UnifiedDataLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top) {

}

template <typename Dtype>
void UnifiedDataLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top){
  for (int i = 0; i < PREFETCH_COUNT; i++) {
    this->prefetch_[i].data_.resize(top.size() - 1);
    for (size_t j = 0; j < this->prefetch_[i].data_.size(); ++j) this->prefetch_[i].data_[j].reset(new Blob<Dtype>);
  }
  DataLayerSetUp(bottom, top);

  // Now, start the prefetch thread. Before calling prefetch, we make two
  // cpu_data calls so that the prefetch thread does not accidentally make
  // simultaneous cudaMalloc calls when the main thread is running. In some
  // GPUs this seems to cause failures if we do not so.
  for (int i = 0; i < PREFETCH_COUNT; i++) {
    for (size_t j = 0; j < prefetch_[i].data_.size() - 1; ++j) {
      if (prefetch_[i].data_[j]->SparseBlob().get() != nullptr)
        prefetch_[i].data_[j]->SparseBlob()->mutable_cpu_data();
      else
        prefetch_[i].data_[j]->mutable_cpu_data();
    }
    prefetch_[i].label_.mutable_cpu_data();
#ifndef CPU_ONLY
      if (Caffe::mode() == Caffe::GPU) {
        for (size_t j = 0; j < prefetch_[i].data_.size() - 1; ++j) {
          if (prefetch_[i].data_[j]->SparseBlob().get() == nullptr)
            prefetch_[i].data_[j]->mutable_gpu_data();
        }
        prefetch_[i].label_.mutable_gpu_data();
      }
#endif

  }
  DLOG(INFO) << "Initializing prefetch";
  StartInternalThread();
  DLOG(INFO) << "Prefetch initialized.";
}

template<class Dtype>
string GetDataDimDesc(const Blob<Dtype>& blob){
  stringstream sstream;
  SparseDataDimension data_dim;
  if (blob.SparseBlob().get() == nullptr){
    sstream << blob.num() << ","
      << blob.channels() << "," << blob.height() << ","
      << blob.width();
  }
  else{
    data_dim = blob.SparseBlob()->dim();
    sstream << "sparse: "
      << "batch: " << data_dim.batch_size << ","
      << "#words: " << data_dim.num_words << ","
      << "#fea: " << data_dim.num_features << ","
      << "vocab: " << data_dim.vocab_size;

  }
  return sstream.str();
}

template <typename Dtype>
void UnifiedDataLayer<Dtype>::DataLayerSetUp(const vector<Blob<Dtype>*>& bottom,
  const vector<Blob<Dtype>*>& top) {
  param_ = this->layer_param_.unified_data_param();
  batch_size_ = this->layer_param_.unified_data_param().batch_size();
  num_output_data_ = top.size() - 1;

  CHECK(this->layer_param_.unified_data_param().data_desc().size() == num_output_data_) << "output blob count mismatch, #outputs: " << num_output_data_ << ", #data desc: " << this->layer_param_.unified_data_param().data_desc().size();

  // Read a data point, and use it to initialize the top blob.
  UnifiedSample& sample = *(reader_.full().pop());

  // Get loader, reshape top and prefect data space
  data_loader_vec_.resize(num_output_data_);
  for (size_t i = 0; i < num_output_data_; ++i){
    data_loader_vec_[i].Init(this->layer_param_.unified_data_param().data_desc(i), phase_);
    data_loader_vec_[i].Reshape(batch_size_, sample.data[i], *top[i]);
    for (size_t prefetch_idx = 0; prefetch_idx < PREFETCH_COUNT; prefetch_idx++) {
      data_loader_vec_[i].Reshape(batch_size_, sample.data[i], *(prefetch_[prefetch_idx].data_[i]));
    }
    LOG(INFO) << "output data size, " << i << ": " << GetDataDimDesc(*top[i]);
  }
  vector<int> label_shape{ batch_size_ };
  top[top.size() - 1]->Reshape(label_shape);
  for (size_t prefetch_idx = 0; prefetch_idx < PREFETCH_COUNT; prefetch_idx++)
    prefetch_[prefetch_idx].label_.Reshape(label_shape);
  
  // push back to free
  reader_.free().push(const_cast<UnifiedSample*>(&sample));

  batch_samples_.resize(num_output_data_);
  batch_images_.resize(num_output_data_);
}

template <typename Dtype>
void CopyBlob(const Blob<Dtype>& src, Blob<Dtype>& dst){
  if (src.SparseBlob().get() == nullptr){
    dst.ReshapeLike(src);
    // Copy the data
    caffe_copy(src.count(), src.cpu_data(),
      dst.mutable_cpu_data());
  }
  else{
    if (dst.SparseBlob().get() == nullptr){
      dst.SparseBlob().reset(new SparseBlob<Dtype>);
    }
    dst.SparseBlob()->Reshape(src.SparseBlob()->dim());
    caffe_copy(src.SparseBlob()->count(), src.SparseBlob()->cpu_data(),
      dst.SparseBlob()->mutable_cpu_data());
  }
}

template <typename Dtype>
void UnifiedDataLayer<Dtype>::Forward_cpu(
  const vector<Blob<Dtype>*>& bottom, const vector<Blob<Dtype>*>& top) {
  UnifiedBatch<Dtype>* batch = prefetch_full_.pop("Data layer prefetch queue empty");
  // Reshape to loaded data.
  for (size_t i = 0; i < top.size() - 1; ++i){
    // Copy the data
    CopyBlob(*batch->data_[i], *top[i]);
  }
  DLOG(INFO) << "Prefetch copied";
  // Copy the labels.
  caffe_copy(batch->label_.count(), batch->label_.cpu_data(),
    top[top.size() - 1]->mutable_cpu_data());

  prefetch_free_.push(batch);
}

template <typename Dtype>
void UnifiedDataLayer<Dtype>::InternalThreadEntry() {
  for (size_t i = 0; i < batch_samples_.size(); ++i){
    batch_samples_[i].resize(batch_size_);
    batch_images_[i].resize(batch_size_);
  }
#ifndef CPU_ONLY
    cudaStream_t stream;
    if (Caffe::mode() == Caffe::GPU) {
      CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    }
#endif

  try {
    while (!must_stop()) {
      UnifiedBatch<Dtype>* batch = prefetch_free_.pop();
      load_batch(batch);
#ifndef CPU_ONLY
        // Why do that??
        if (Caffe::mode() == Caffe::GPU) {
          // sync each blob
          for (int i = 0; i < batch->data_.size(); i++) {
            if (batch->data_[i]->SparseBlob().get() == nullptr) {
              batch->data_[i]->data().get()->async_gpu_push(stream);
              CUDA_CHECK(cudaStreamSynchronize(stream));
            }
          }
        }
#endif
        prefetch_full_.push(batch);
      }
    }
    catch (boost::thread_interrupted&) {
      // Interrupted exception is expected on shutdown
    }
#ifndef CPU_ONLY
    if (Caffe::mode() == Caffe::GPU) {
      CUDA_CHECK(cudaStreamDestroy(stream));
    }
#endif
}

// This function is called on prefetch thread
//  If input contains file path to IUB files, this function will also copy the prefetched images
template<typename Dtype>
void UnifiedDataLayer<Dtype>::load_batch(UnifiedBatch<Dtype>* batch) {
  CPUTimer batch_timer;
  batch_timer.Start();
  double read_time = 0;
  double trans_time = 0;
  CPUTimer timer;
  for (int i = 0; i < num_output_data_; i++) {
    if (batch->data_[i]->SparseBlob().get() == nullptr) {
      CHECK(batch->data_[i]->count()) << "blob count should not be zero.";
    }
    else {
      CHECK(batch->data_[i]->SparseBlob()->count()) << "sparse blob count should not be zero.";
    }
  }

  // Reshape according to the first datum of each batch
  // on single input batches allows for inputs of varying dimension.
  const int batch_size = this->layer_param_.data_param().batch_size();
  Dtype* prefetch_label = batch->label_.mutable_cpu_data();
  // copy from reader to batch samples, then loading
  for (int num = 0; num < batch_size_; ++num){
    timer.Start();

    // get a sample
	UnifiedSample& sample = *(reader_.full().pop("Waiting for data"));
    read_time += timer.MicroSeconds();

    // Copy to batch_samples_
    timer.Start();
    for (size_t i = 0; i < batch_samples_.size(); ++i){
      batch_samples_[i][num] = sample.data[i];
      // Copy cached image from UnifiedSample into batch_images_
      if (this->layer_param_.unified_data_param().data_desc(i).type() == DataDescriptionParameter::IMAGEIUB)
        batch_images_[i][num] = (*(sample.image_data))[i];
    }

    // Copy Label
    prefetch_label[num] = (Dtype)(sample.label);
    trans_time += timer.MicroSeconds();
	reader_.free().push(const_cast<UnifiedSample*>(&sample));
  }

  // Loading
  for (size_t i = 0; i < batch_samples_.size(); ++i){
    if (this->layer_param_.unified_data_param().data_desc(i).type() == DataDescriptionParameter::IMAGEIUB)
      data_loader_vec_[i].Load(batch_images_[i], *batch->data_[i]); // Only transformation is done 
    else
      data_loader_vec_[i].Load(batch_samples_[i], *batch->data_[i]);
  }

  timer.Stop();
  batch_timer.Stop();
  DLOG(INFO) << "Prefetch batch: " << batch_timer.MilliSeconds() << " ms.";
  DLOG(INFO) << "     Read time: " << read_time / 1000 << " ms.";
  DLOG(INFO) << "Transform time: " << trans_time / 1000 << " ms.";
}

INSTANTIATE_CLASS(UnifiedDataLayer);
REGISTER_LAYER_CLASS(UnifiedData);
}