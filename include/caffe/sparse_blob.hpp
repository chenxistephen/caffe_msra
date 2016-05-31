#ifndef CAFFE_SPARSE_BLOB_HPP_
#define CAFFE_SPARSE_BLOB_HPP_

#include <algorithm>
#include <string>
#include <vector>

#include "caffe/common.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/syncedmem.hpp"
#include "caffe/util/math_functions.hpp"

namespace caffe 
{
	/**
	* @brief 
	*/
	template <typename Dtype>
    struct SparseItem
    {
        static const uint32_t NULL_ID = uint32_t(-1);

        uint32_t sample_id;
        uint32_t word_id;
        uint32_t feature_id;
        Dtype val;

        SparseItem(uint32_t sample_id_, uint32_t word_id_, uint32_t feature_id_, Dtype val_)
            :sample_id(sample_id_), word_id(word_id_), feature_id(feature_id_), val(val_) {

        }

        SparseItem(){
            sample_id = NULL_ID; word_id = NULL_ID; feature_id = NULL_ID; val = 0;
        }

        bool IsNull() const {
            return sample_id == NULL_ID || word_id == NULL_ID || feature_id == NULL_ID;
        }

        static SparseItem NullItem() {
            static SparseItem item;
            item.sample_id = NULL_ID; item.word_id = NULL_ID; item.feature_id = NULL_ID; item.val = 0;
            return item;
        }
	};


    struct SparseDataDimension{
        uint32_t batch_size;
        uint32_t num_words;
        uint32_t num_features;
        uint32_t vocab_size;

        SparseDataDimension()
            :batch_size(0), num_words(0), num_features(0), vocab_size(0) {}

        SparseDataDimension(uint32_t batch_size_, uint32_t num_words_, uint32_t num_features_, uint32_t vocab_size_)
            :batch_size(batch_size_), num_words(num_words_), num_features(num_features_), vocab_size(vocab_size_) {}
    };

	/**
	* @brief
	*/
	template <typename Dtype>
	class SparseBlob {
	public:
		typedef SparseItem<Dtype> t_item;

		SparseBlob() : data_(), diff_(), count_(0), capacity_(0) {}

		/// @brief In text embedding, num is the number of text documents in a batch, height is the maximum number of words in a document, width is the maximum number of non-zero elements in the quantized vector of a word. 
		///        If id=-1 then it means the end of this sparse vector
		///        Here the word means the basic unit to input to Caffe for training (e.g., in DSSM it is a English word, in "Text Understanding from Scratch" it is a letter)
        explicit SparseBlob(const int batch_size, const int num_words, const int num_features, const int vocab_size);
        void Reshape(const int batch_size, const int num_words, const int num_features, const int vocab_size);
        void Reshape(const SparseDataDimension& dim);

		const t_item* cpu_data() const;
		void set_cpu_data(t_item* data);
		const t_item* gpu_data() const;
		const t_item* cpu_diff() const;
		const t_item* gpu_diff() const;
		t_item* mutable_cpu_data();
		t_item* mutable_gpu_data();
		t_item* mutable_cpu_diff();
		t_item* mutable_gpu_diff();
        inline int count() const { return count_; }
        inline int vocab_size() const { return vocab_size_; }
        inline int word_count() const { return shape_[1]; }
        inline int batch_size() const { return shape_[0]; }
        inline SparseDataDimension dim() const { return SparseDataDimension(batch_size(), word_count(), shape_[2], vocab_size()); }

		//void Update();
		//void FromProto(const SparseBlobProto& proto, bool reshape = true);
		//void ToProto(SparseBlobProto* proto, bool write_diff = false) const;
    protected:
        shared_ptr<SyncedMemory> data_;
        shared_ptr<SyncedMemory> diff_;
        vector<int> shape_;
        int count_;
        int capacity_;
        int vocab_size_;

        DISABLE_COPY_AND_ASSIGN(SparseBlob);
	};

}  // namespace caffe

#endif  // CAFFE_SPARSE_BLOB_HPP_
