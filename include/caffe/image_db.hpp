#pragma once

#include <string>
#include <utility>
#include <vector>

#include "caffe/proto/caffe.pb.h"
#include "caffe/blob.hpp"
#include "caffe/common.hpp"

namespace caffe {
    struct Rect {
        float left = 0.f;
        float top = 0.f;
        float width = 0.f;
        float height = 0.f;

        bool IsNull() const {
            return width == 0.f || height == 0.f;
        }
    };

    struct ImageData {
        vector<uint8_t> data;
        Rect bbox;
        std::string key = "";
    };

    class IImageDB {
    public:
        virtual void Init(const ImageDBParameter& param) = 0;
        virtual void Read(const string& key, ImageData& data) = 0;

    public:
		static shared_ptr<IImageDB> CreateInstanceAndInit(const DataDescriptionParameter& param);
    };

    class FileSystemImageDB : public IImageDB {
    public:
        void Init(const ImageDBParameter& param);
        virtual void Read(const string& key, ImageData& data);
    protected:
        void GetImageInfo(const string& key, string& imgage_key, Rect& crop);
    private:
        string root_folder_;
    };

	class IUBFileSystemImageDB : public IImageDB {
	public:
		void Init(const ImageDBParameter& param);
		virtual void Read(const string& key, ImageData& data);
	protected:
		void GetImageInfo(const string& key, string& imgage_key, Rect& crop);
	private:
		std::string root_folder_;
		std::string current_filename;
		std::ifstream src;
	};

}