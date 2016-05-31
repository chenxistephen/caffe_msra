#include "caffe/data_loader.hpp"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "caffe/util/benchmark.hpp"
#include <boost/timer/timer.hpp>
#include "ppl.h"
using namespace concurrency;
#include "caffe/base64.hpp"

namespace caffe {

    template<class Dtype>
    void TextDataLoader<Dtype>::Init(const DataDescriptionParameter& param, Phase phase){
        param_ = param;
        phase_ = phase;
        CHECK(param_.type() == DataDescriptionParameter::TEXT) << "data loader type mismatch";
        text_hashing_algo_ = ITextHashing<Dtype>::CreateInstanceAndInit(param_.text_hashing_param().algo_name(), param_.text_hashing_param());
    }


    template<class Dtype>
    void TextDataLoader<Dtype>::Load(const vector<string>& input_data, Blob<Dtype>& out_blob){
        if (out_blob.SparseBlob().get() == nullptr) {
            out_blob.SparseBlob().reset(new SparseBlob<Dtype>);
        }
        text_hashing_algo_->Extract(input_data, *out_blob.SparseBlob());
    }

    template<class Dtype>
    void TextDataLoader<Dtype>::Reshape(int batch_size, string input_sample, Blob<Dtype>& out_blob){
        if (out_blob.SparseBlob().get() == nullptr) {
            out_blob.SparseBlob().reset(new SparseBlob<Dtype>);
        }
        out_blob.SparseBlob()->Reshape(text_hashing_algo_->GetDataDimension(batch_size));
    }

    
    template<class Dtype>
    void ImageDataLoader<Dtype>::Init(const DataDescriptionParameter& param, Phase phase){
        param_ = param;
        phase_ = phase;
        CHECK(param_.type() == DataDescriptionParameter::IMAGE) << "data loader type mismatch";
        image_db_ = IImageDB::CreateInstanceAndInit(param_);
        image_transformer_ = IImageDataTransformer<Dtype>::CreateInstanceAndInit(param_.image_data_desc_param().transform_param().type(), param_.image_data_desc_param().transform_param(), phase);
    }


    template<class Dtype>
    void ImageDataLoader<Dtype>::Load(const vector<string>& input_data, Blob<Dtype>& out_blob){
        vector<ImageData> image_data(input_data.size());
        CPUTimer batch_timer;
        image_data.resize(input_data.size());
        batch_timer.Start();
		//for (size_t i = 0; i < input_data.size(); ++i) {
		parallel_for((size_t)0, input_data.size(), [&](size_t i){
            image_db_->Read(input_data[i], image_data[i]);
        });
        batch_timer.Stop();
        DLOG(INFO) << "image db load: " << batch_timer.MilliSeconds() << " ms.";

        batch_timer.Start();
		image_transformer_->Transform(image_data, out_blob);
        batch_timer.Stop();
        DLOG(INFO) << "image data transform: " << batch_timer.MilliSeconds() << " ms.";
    }

    template<class Dtype>
    void ImageDataLoader<Dtype>::Reshape(int batch_size, string input_sample, Blob<Dtype>& out_blob){
        Blob<Dtype> tmp_blob;
        vector<int> dim;
        Load(vector < string > {input_sample}, tmp_blob);
        dim = tmp_blob.shape();
        dim[0] = batch_size;
        out_blob.Reshape(dim);
    }

    template<class Dtype>
    void FeatureDataLoader<Dtype>::Init(const DataDescriptionParameter& param, Phase phase){
        param_ = param;
        phase_ = phase;
        CHECK(param_.type() == DataDescriptionParameter::FEATURE) << "data loader type mismatch";
    }

    template<class Dtype>
    void FeatureDataLoader<Dtype>::Load(const vector<string>& input_data, Blob<Dtype>& out_blob){
        int dim = 0;
        vector<int> idx(2);
        for (size_t i = 0; i < input_data.size(); ++i){
            vector<string> ps;
            string feature_string = input_data[i];
            boost::trim(feature_string);
            boost::split(ps, feature_string, boost::is_any_of(" "));
            if (i == 0){
                out_blob.Reshape(vector < int > {(int)input_data.size(), (int)ps.size()});
                dim = ps.size();
            }
            CHECK(ps.size() == dim) << "feature dimension mismatch, prev: " << dim << "; cur: " << ps.size();
            idx[0] = (int)i; idx[1] = 0;
            Dtype* p_data = out_blob.mutable_cpu_data() + out_blob.offset(idx);
            
            for (size_t j = 0; j < ps.size(); ++j){
                try {
                  p_data[j] = std::stof(ps[j]);
                }
                catch (boost::bad_lexical_cast const&) {
                    CHECK(false) << "the feature should be float, part: " << ps[j] << ", line: " << input_data[i];
                }
            }
        }
    }

    template<class Dtype>
    void FeatureDataLoader<Dtype>::Reshape(int batch_size, string input_sample, Blob<Dtype>& out_blob){
        Blob<Dtype> tmp_blob;
        vector<int> dim;
        Load(vector < string > {input_sample}, tmp_blob);
        dim = tmp_blob.shape();
        dim[0] = batch_size;
        out_blob.Reshape(dim);
    }

    template<class Dtype>
    void ImageBase64DataLoader<Dtype>::Init(const DataDescriptionParameter& param, Phase phase){
        param_ = param;
        phase_ = phase;
        CHECK(param_.type() == DataDescriptionParameter::IMAGEBASE64) << "data loader type mismatch";
        image_transformer_ = IImageDataTransformer<Dtype>::CreateInstanceAndInit(param_.image_data_desc_param().transform_param().type(), param_.image_data_desc_param().transform_param(), phase);
    }

    template<class Dtype>
    void ImageBase64DataLoader<Dtype>::Load(const vector<string>& input_data, Blob<Dtype>& out_blob){
        CPUTimer batch_timer;
        vector<ImageData> image_data(input_data.size());
        image_data.resize(input_data.size());

        batch_timer.Start();
        parallel_for((size_t)0, input_data.size(), [&](size_t i){
            IU::Base64Decode(input_data[i], image_data[i].data);
        }
        );
        batch_timer.Stop();
        DLOG(INFO) << "image data base64 decoding: " << batch_timer.MilliSeconds() << " ms.";
        batch_timer.Start();
        //for (size_t i = 0; i < input_data.size(); ++i){
        image_transformer_->Transform(image_data, out_blob);
        batch_timer.Stop();
        DLOG(INFO) << "image data transform: " << batch_timer.MilliSeconds() << " ms.";
    }

    template<class Dtype>
    void ImageBase64DataLoader<Dtype>::Reshape(int batch_size, string input_sample, Blob<Dtype>& out_blob){
        Blob<Dtype> tmp_blob;
        vector<int> dim;
        Load(vector < string > {input_sample}, tmp_blob);
        dim = tmp_blob.shape();
        dim[0] = batch_size;
        out_blob.Reshape(dim);
    }

	/**
	* IUB Data loader
	*/
	template<class Dtype>
	void ImageIUBDataLoader<Dtype>::Init(const DataDescriptionParameter& param, Phase phase){
		param_ = param;
		phase_ = phase;
		CHECK(param_.type() == DataDescriptionParameter::IMAGEIUB) << "data loader type mismatch";
		image_db_ = IImageDB::CreateInstanceAndInit(param_);
		image_transformer_ = IImageDataTransformer<Dtype>::CreateInstanceAndInit(param_.image_data_desc_param().transform_param().type(), param_.image_data_desc_param().transform_param(), phase);
	}

	template<class Dtype>
	void ImageIUBDataLoader<Dtype>::Load(const vector<string>& input_data, Blob<Dtype>& out_blob){
		vector<ImageData> image_data(input_data.size());
		CPUTimer batch_timer;

		image_data.resize(input_data.size());

		batch_timer.Start();
		for (size_t i = 0; i < input_data.size(); ++i) {
			image_db_->Read(input_data[i], image_data[i]);
		}
		batch_timer.Stop();
		DLOG(INFO) << "image db load: " << batch_timer.MilliSeconds() << " ms.";

		batch_timer.Start();
		image_transformer_->Transform(image_data, out_blob);
		batch_timer.Stop();
		DLOG(INFO) << "image data transform: " << batch_timer.MilliSeconds() << " ms.";
	}

	template<class Dtype>
	void ImageIUBDataLoader<Dtype>::Load(const vector<ImageData>& input_image_data, Blob<Dtype>& out_blob){
		CPUTimer batch_timer;
		batch_timer.Start();
		image_transformer_->Transform(input_image_data, out_blob);
		batch_timer.Stop();
		DLOG(INFO) << "image data transform: " << batch_timer.MilliSeconds() << " ms.";
	}

	template<class Dtype>
	void ImageIUBDataLoader<Dtype>::Reshape(int batch_size, string input_sample, Blob<Dtype>& out_blob){
		Blob<Dtype> tmp_blob;
		vector<int> dim;
		Load(vector < string > {input_sample}, tmp_blob);
		dim = tmp_blob.shape();
		dim[0] = batch_size;
		out_blob.Reshape(dim);
	}


    template<class Dtype>
    void UnifiedDataLoader<Dtype>::Init(const DataDescriptionParameter& param, Phase phase){
        switch (param.type()) {
        case DataDescriptionParameter::TEXT:
            m_data_loader.reset(new TextDataLoader<Dtype>);
            break;
        case DataDescriptionParameter::IMAGE:
            m_data_loader.reset(new ImageDataLoader<Dtype>);
            break;
        case DataDescriptionParameter::FEATURE:
            m_data_loader.reset(new FeatureDataLoader<Dtype>);
            break;
        case DataDescriptionParameter::IMAGEBASE64:
            m_data_loader.reset(new ImageBase64DataLoader<Dtype>);
            break;
		case DataDescriptionParameter::IMAGEIUB:
			m_data_loader.reset(new ImageIUBDataLoader<Dtype>);
			break;
        default:
            LOG(FATAL) << "non-supported type";
        }
        m_data_loader->Init(param, phase);
    }

    template<class Dtype>
    void UnifiedDataLoader<Dtype>::Load(const vector<string>& input_data, Blob<Dtype>& out_blob){
        CHECK(m_data_loader.get() != nullptr) << "data loader not initialized";
        m_data_loader->Load(input_data, out_blob);
    }

	template<class Dtype>
	void UnifiedDataLoader<Dtype>::Load(const vector<ImageData>& input_image_data, Blob<Dtype>& out_blob){
		CHECK(m_data_loader.get() != nullptr) << "data loader not initialized";
		dynamic_cast<ImageIUBDataLoader<Dtype>*>(m_data_loader.get())->Load(input_image_data, out_blob);
		//dynamic_pointer_cast<shared_ptr<ImageIUBDataLoader<Dtype> > >(m_data_loader)->Load(input_image_data, out_blob);
	}

    template<class Dtype>
    void UnifiedDataLoader<Dtype>::Reshape(int batch_size, string input_sample, Blob<Dtype>& out_blob){
        CHECK(m_data_loader.get() != nullptr) << "data loader not initialized";
        m_data_loader->Reshape(batch_size, input_sample, out_blob);
    }

    INSTANTIATE_CLASS(UnifiedDataLoader);
}