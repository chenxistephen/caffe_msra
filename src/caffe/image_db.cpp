#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

#include "caffe/image_db.hpp"

using namespace std;
using namespace boost::algorithm;

namespace caffe {

	shared_ptr<IImageDB> IImageDB::CreateInstanceAndInit(const DataDescriptionParameter& param) {
        shared_ptr<IImageDB> obj;
		if (param.image_data_desc_param().image_db_param().type() == "filesystem" && param.type() == DataDescriptionParameter::IMAGE)
			obj.reset(new FileSystemImageDB());
		else if (param.image_data_desc_param().image_db_param().type() == "IUBfilesystem" && param.type() == DataDescriptionParameter::IMAGEIUB)
			obj.reset(new IUBFileSystemImageDB());
		CHECK(obj.get() != nullptr) << "non-supported algorithm: " << param.image_data_desc_param().image_db_param().type();
		obj->Init(param.image_data_desc_param().image_db_param());
        return obj;
    }

    void FileSystemImageDB::Init(const ImageDBParameter& param) {
        root_folder_ = param.source();
        if (root_folder_.length() > 0 && root_folder_[root_folder_.length() - 1] != '\\' && root_folder_[root_folder_.length() - 1] != '/')
            root_folder_.append("\\");
    }

    void FileSystemImageDB::GetImageInfo(const string& key, string& imgage_key, Rect& bbox) {
        vector <string> ps;
        // Use "|" to split bbox and file name. We can't use ":" for absolute path may contain ":"
        //boost::split(ps, key, boost::is_any_of(":"));
        boost::split(ps, key, boost::is_any_of("|"));
        imgage_key = ps[0];
        if (ps.size() >= 2) {
            vector<string> bbps;
            boost::split(bbps, ps[1], boost::is_any_of(", "));
            CHECK_EQ(bbps.size(), 4) << "there must be 4 numbers for bbox, now " << ps[1];
            try {
                bbox.left = boost::lexical_cast<float>(bbps[0]);
                bbox.top = boost::lexical_cast<float>(bbps[1]);
                bbox.width = boost::lexical_cast<float>(bbps[2]);
                bbox.height = boost::lexical_cast<float>(bbps[3]);
            }
            catch (boost::bad_lexical_cast const&) {
                LOG(ERROR) << "incorrect line, bbox should be all numbers!. It is now: " << ps[1];
            }
        }
    }

    void ReadFileIntoBuffer(const char *filename, vector<uint8_t>& data)
    {
#ifndef NDEBUG
        //LOG(INFO) << "processing " << filename;
#endif
        std::ifstream file(filename, std::ios::binary);
        CHECK(file.good()) << "file can't be read: " << filename << "!!!!! " << file.bad() << " " << file.eof() << " " << file.fail();

        file.unsetf(std::ios::skipws);

        file.seekg(0, std::ios::end);
        std::streampos file_size = file.tellg();
        file.seekg(0, std::ios::beg);

        data.reserve(file_size);
        data.insert(data.begin(),
            std::istream_iterator<uint8_t>(file),
            std::istream_iterator<uint8_t>());
    }

    void FileSystemImageDB::Read(const string& key, ImageData& data) {
        string image_key;
        GetImageInfo(key, image_key, data.bbox);
        data.key = key;
        string filename = root_folder_ + image_key;
        ReadFileIntoBuffer(filename.c_str(), data.data);
    }

	/*
		IUBFilesystem
	*/
	void IUBFileSystemImageDB::Init(const ImageDBParameter& param) {
		root_folder_ = param.source();
		if (root_folder_.length() > 0 && root_folder_[root_folder_.length() - 1] != '\\' && root_folder_[root_folder_.length() - 1] != '/')
			root_folder_.append("\\");
	}

	void IUBFileSystemImageDB::GetImageInfo(const string& key, string& imgage_key, Rect& bbox) {
		vector <string> ps;
		// Use "|" to split bbox and file name. We can't use ":" for absolute path may contain ":"
		//boost::split(ps, key, boost::is_any_of(":"));
		boost::split(ps, key, boost::is_any_of("|"));
		imgage_key = ps[0];
		if (ps.size() >= 2) {
			vector<string> bbps;
			boost::split(bbps, ps[1], boost::is_any_of(", "));
			CHECK_EQ(bbps.size(), 4) << "there must be 4 numbers for bbox, now " << ps[1];
			try {
				bbox.left = boost::lexical_cast<float>(bbps[0]);
				bbox.top = boost::lexical_cast<float>(bbps[1]);
				bbox.width = boost::lexical_cast<float>(bbps[2]);
				bbox.height = boost::lexical_cast<float>(bbps[3]);
			}
			catch (boost::bad_lexical_cast const&) {
				LOG(ERROR) << "incorrect line, bbox should be all numbers!. It is now: " << ps[1];
			}
		}
	}

	void IUBFileSystemImageDB::Read(const string& key, ImageData& data) {
		string image_key;
		GetImageInfo(key, image_key, data.bbox);
		data.key = key;
		vector<string> image_file_tuple;
		boost::split(image_file_tuple, image_key, boost::is_any_of(";"));
		string filename = image_file_tuple[0];
		unsigned long long offset = std::stoull(image_file_tuple[1]);
		unsigned long long size = std::stoull(image_file_tuple[2]);
		string fullpath = root_folder_ + filename;
		if (current_filename.empty() || current_filename.compare(filename) != 0) {
			src.close();
			current_filename = fullpath;
			src.open(current_filename.c_str(), ios_base::binary);
		}

		data.data.resize(size);
		if (src.cur != src.beg + offset)
			src.seekg(offset, src.beg);
		src.read((char*)&data.data[0], size);
	}
}