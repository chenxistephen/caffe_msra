#include "caffe/sparse_blob.hpp"
#include "caffe/blob.hpp"

namespace caffe{

template <typename Dtype>
void SparseBlob<Dtype>::Reshape(const int batch_size, const int num_words, const int num_features, const int vocab_size) {
    vocab_size_ = vocab_size;
	shape_.resize(3);
    shape_[0] = batch_size;
    shape_[1] = num_words;
    shape_[2] = num_features;
    count_ = batch_size * num_words * num_features;
	if (count_ > capacity_) {
		capacity_ = count_;
		data_.reset(new SyncedMemory(capacity_ * sizeof(SparseItem<Dtype>)));
		diff_.reset(new SyncedMemory(capacity_ * sizeof(SparseItem<Dtype>)));
	}
}

template <typename Dtype>
void SparseBlob<Dtype>::Reshape(const SparseDataDimension& dim) {
    Reshape((int)dim.batch_size, (int)dim.num_words, (int)dim.num_features, (int)dim.vocab_size);
}

template <typename Dtype>
SparseBlob<Dtype>::SparseBlob(const int batch_size, const int num_words, const int num_features, const int vocab_size)
	// capacity_ must be initialized before calling Reshape
	: capacity_(0) {
    Reshape(batch_size, num_words, num_features, vocab_size);
}

template <typename Dtype>
const SparseItem<Dtype>* SparseBlob<Dtype>::cpu_data() const {
	CHECK(data_);
	return (const SparseItem<Dtype>*)data_->cpu_data();
}

template <typename Dtype>
void SparseBlob<Dtype>::set_cpu_data(SparseItem<Dtype>* data) {
	CHECK(data);
	data_->set_cpu_data(data);
}

template <typename Dtype>
const SparseItem<Dtype>* SparseBlob<Dtype>::gpu_data() const {
	CHECK(data_);
	return (const SparseItem<Dtype>*)data_->gpu_data();
}

template <typename Dtype>
const SparseItem<Dtype>* SparseBlob<Dtype>::cpu_diff() const {
	CHECK(diff_);
	return (const SparseItem<Dtype>*)diff_->cpu_data();
}

template <typename Dtype>
const SparseItem<Dtype>* SparseBlob<Dtype>::gpu_diff() const {
	CHECK(diff_);
	return (const SparseItem<Dtype>*)diff_->gpu_data();
}

template <typename Dtype>
SparseItem<Dtype>* SparseBlob<Dtype>::mutable_cpu_data() {
	CHECK(data_);
	return static_cast<SparseItem<Dtype>*>(data_->mutable_cpu_data());
}

template <typename Dtype>
SparseItem<Dtype>* SparseBlob<Dtype>::mutable_gpu_data() {
	CHECK(data_);
	return static_cast<SparseItem<Dtype>*>(data_->mutable_gpu_data());
}

template <typename Dtype>
SparseItem<Dtype>* SparseBlob<Dtype>::mutable_cpu_diff() {
	CHECK(diff_);
	return static_cast<SparseItem<Dtype>*>(diff_->mutable_cpu_data());
}

template <typename Dtype>
SparseItem<Dtype>* SparseBlob<Dtype>::mutable_gpu_diff() {
	CHECK(diff_);
	return static_cast<SparseItem<Dtype>*>(diff_->mutable_gpu_data());
}

template<class Dtype>
void SparseBlob2DenseBlob(const SparseBlob<Dtype>& input_sprase_blob, Blob<Dtype>& dense_blob){
    const SparseItem<Dtype>* p_sprase_items = input_sprase_blob.cpu_data();
    dense_blob.Reshape(input_sprase_blob.batch_size(), input_sprase_blob.vocab_size(), input_sprase_blob.word_count(), 1);
    Dtype* p_dense = dense_blob.mutable_cpu_data();
    caffe_set(dense_blob.count(), Dtype(0.0), p_dense);
    for (size_t i = 0; i < input_sprase_blob.count(); ++i){
        if (p_sprase_items[i].IsNull()) break;
        p_dense[dense_blob.offset(p_sprase_items[i].sample_id, p_sprase_items[i].feature_id, p_sprase_items[i].word_id, 0)] = p_sprase_items[i].val;
    }
}

template<class Dtype>
void DenseBlob2SparseBlob(const Blob<Dtype>& input_dense_blob, SparseBlob<Dtype>& sparse_blob){
    CHECK((input_dense_blob.num_axes() == 4 && input_dense_blob.shape()[3] == 1))<< "input must have 4 axes and the width is 1.";
    auto shape = input_dense_blob.shape();
    int num_words = shape[2];
    int vocab_size = shape[1];
    sparse_blob.Reshape(shape[0], num_words, vocab_size, vocab_size);
    SparseItem<Dtype>* p_sprase_items = sparse_blob.mutable_cpu_data();
    const Dtype* p_dense = input_dense_blob.cpu_data();
    for (int i = 0; i < sparse_blob.count(); ++i){
        p_sprase_items[i] = SparseItem<Dtype>::NullItem();
    }
    int m = 0;
    for (int i = 0; i < shape[0]; ++i){
        for (int j = 0; j < shape[1]; ++j){
            for (int k = 0; k < shape[2]; ++k){
                auto val = p_dense[input_dense_blob.offset(i, j, k, 0)];
                if (val != Dtype(0)){
                    p_sprase_items[m].sample_id = i;
                    p_sprase_items[m].feature_id = j;
                    p_sprase_items[m].word_id = k;
                    p_sprase_items[m].val = val;
                        ++m;
                }
            }
        }
    }
}

template void SparseBlob2DenseBlob<float>(const SparseBlob<float>& input_sprase_blob, Blob<float>& dense_blob);
template void DenseBlob2SparseBlob<float>(const Blob<float>& input_dense_blob, SparseBlob<float>& sprase_blob);
template void SparseBlob2DenseBlob<double>(const SparseBlob<double>& input_sprase_blob, Blob<double>& dense_blob);
template void DenseBlob2SparseBlob<double>(const Blob<double>& input_dense_blob, SparseBlob<double>& sprase_blob);

INSTANTIATE_CLASS(SparseBlob);
template class SparseBlob<int>;
template class SparseBlob<unsigned int>;

}