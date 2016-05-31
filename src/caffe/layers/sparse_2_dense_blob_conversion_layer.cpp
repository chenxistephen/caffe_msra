#include "caffe/layers/nlp_layers.hpp"
#include <fstream>

namespace caffe{

    template<typename Dtype>
    void Sparse2DenseBlobConversionLayer<Dtype>::LayerSetUp(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        CHECK(bottom.size() >= 1 && bottom[0]->SparseBlob().get() != nullptr)
            << "The input should be sparse blob";
    }

    template<typename Dtype>
    void Sparse2DenseBlobConversionLayer<Dtype>::Reshape(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        CHECK(top.size() >= 1)
            << "There is at least one output from Sparse2DenseBlobConversionLayer.";
        CHECK(bottom.size() >= 1)
            << "There is at least one input from Sparse2DenseBlobConversionLayer.";
        for (int i = 0; i < bottom.size(); ++i) {
            int count_images = bottom[i]->SparseBlob()->batch_size();
            auto word_count = bottom[i]->SparseBlob()->word_count();
            auto vocab_size = bottom[i]->SparseBlob()->vocab_size();
            vector<int> top_shape{ count_images, vocab_size, word_count, 1 };
            top[i]->Reshape(top_shape);
        }
    }

    template<typename Dtype>
    void WriteBlob2Sparse(const char* filename, const Blob<Dtype>& blob){
        std::ofstream of(filename);
        auto& shape = blob.shape();
        CHECK(shape.size() == 4 && shape[3] == 1) << "shape must be 4 axes and the width is 1.";
        for (int n = 0; n < shape[0]; ++n){
            for (int h = 0; h < shape[2]; ++h){
                for (int c = 0; c < shape[1]; ++c){
                    auto d = blob.data_at(n, c, h, 0);
                    if (d != 0)
                        of << c << ":" << d<<"#";
                }
            }
            of << std::endl;
        }
    }

    template<typename Dtype>
    void Sparse2DenseBlobConversionLayer<Dtype>::Forward_cpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        for (int i = 0; i < bottom.size(); ++i) {
            CHECK(bottom[i]->SparseBlob().get() != nullptr)
                << "The input should be sparse blob";
            auto sparse_input = bottom[i]->SparseBlob();
            SparseBlob2DenseBlob(*sparse_input, *(top[i])); //dense: num_images, vocab_size, #words_in_doc
            //WriteBlob2Sparse("test_blob_convert.txt", *top[i]);
        }
    }

    template<typename Dtype>
    void Sparse2DenseBlobConversionLayer<Dtype>::Forward_gpu(const vector<Blob<Dtype>*>& bottom,
        const vector<Blob<Dtype>*>& top){
        Forward_cpu(bottom, top);
    }

    template<typename Dtype>
    void Sparse2DenseBlobConversionLayer<Dtype>::Backward_cpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom){
    }

    template<typename Dtype>
    void Sparse2DenseBlobConversionLayer<Dtype>::Backward_gpu(const vector<Blob<Dtype>*>& top,
        const vector<bool>& propagate_down, const vector<Blob<Dtype>*>& bottom){
        Backward_cpu(top, propagate_down, bottom);
    }

    INSTANTIATE_CLASS(Sparse2DenseBlobConversionLayer);
    REGISTER_LAYER_CLASS(Sparse2DenseBlobConversion);
}