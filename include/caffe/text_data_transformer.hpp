#ifndef CAFFE_TEXT_DATA_TRANSFORMER_HPP_
#define CAFFE_TEXT_DATA_TRANSFORMER_HPP_

#include <string>
#include <vector>
#include <unordered_map>
#include <boost/tokenizer.hpp>

#include "caffe/common.hpp"
#include "caffe/proto/caffe.pb.h"
#include "caffe/sparse_blob.hpp"

namespace caffe{
    using std::vector;
    using std::string;
    using std::unordered_map;
    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

    template<class Dtype>
    class ITextHashing{
    public:
        virtual ~ITextHashing() {}
        virtual void Init(const TextHashingParameter& param) = 0;
        virtual void Extract(const vector<string> &data, SparseBlob<Dtype>& out_blob) const = 0;
        virtual SparseDataDimension GetDataDimension(uint32_t batch_size) const = 0;

    public:
        static shared_ptr<ITextHashing<Dtype>> CreateInstanceAndInit(const string& algo_name, const TextHashingParameter& param);
    };

    template<class Dtype>
    class TextHashingBase : public ITextHashing<Dtype>{
    public:
        void Init(const TextHashingParameter& param) override{
            param_ = param;
            if (param_.has_build_dict_param()) {
                std::ifstream infile(param_.build_dict_param().save_path());
                CHECK(!infile) << "dictionary has already been build at the given path, please delete the old dictionary file and rebuild dictionary.";
                infile.close();
                BuildDictionary();
            }
            else {
                LoadDictionary();
            }
        }

    protected:
        TextHashingParameter param_;

    protected:
        virtual void LoadDictionary() = 0;
        virtual void BuildDictionary() = 0;
    };

    template<class Dtype>
    class WordHashingBase : public TextHashingBase<Dtype>{
    public:
        SparseDataDimension GetDataDimension(uint32_t batch_size) const {
            return SparseDataDimension(batch_size, param_.max_num_words_per_doc(), param_.max_num_non_zeros_per_word(), dict_.size());
        }
    protected:
        void LoadDictionary() override;
        void BuildDictionary() override;
        /* !
         * \brief parse the given text to the desired representation units (word/char gram)
         * \param str the string to be parsed
         * \return representation units in a vector
        */
        virtual vector<string> ParseText(string& str) const = 0;
        /* !
        * \brief split the str into n_word_grams using delimiters
        * \param str the string to be tokenized
        * \param ngram for each value i in ngram, i-gram will be extracted
        * \return tokens in a vector
        */
        vector<string>  Tokenize(string& str, vector<uint32_t>& ngram) const;
    protected:
        unordered_map<string, int> dict_;
        boost::char_separator<char> delimiters_{"\x20\x0A\x0C\x09\x0B\x0D\xA0\x21\x22\x23\x26\x27\x28\x29\x2A\x2B\x2C\x2D\x2E\x2F\x40\x5B\x5C\x5D\x5E\x5F\x60\x7B\x7C\x7D\x7E\x7F"};
    };

    template<class Dtype>
    class ChargramHashing : public TextHashingBase<Dtype>{
    public:
        void Init(const TextHashingParameter& param) override;
        void Extract(const vector<string> &data, SparseBlob<Dtype>& out_blob) const override;
        SparseDataDimension GetDataDimension(uint32_t batch_size) const {
            return SparseDataDimension(batch_size, param_.max_num_words_per_doc(), param_.max_num_non_zeros_per_word(), dict_.size());
        }
    protected:
        void LoadDictionary() override;
        void BuildDictionary() override;
    protected:
        unordered_map<char, int> dict_;
    };

    template<class Dtype>
    class BOWHashing : public WordHashingBase<Dtype> {
    public:
        SparseDataDimension GetDataDimension(uint32_t batch_size) const {
            return SparseDataDimension((int)batch_size, 1, param_.max_num_words_per_doc(), dict_.size());
        }
    protected:
        /* !
        * \brief For BOWHashing, ParseText will return word_n_gram
        */
        virtual vector<string> ParseText(string& str) const override;
    public:
        void Init(const TextHashingParameter& param) override;
        void Extract(const vector<string> &data, SparseBlob<Dtype>& out_blob) const override;
    };

    template<class Dtype>
    class NLettergramHashing : public WordHashingBase<Dtype> {
    public:
        SparseDataDimension GetDataDimension(uint32_t batch_size) const {
            return SparseDataDimension((int)batch_size, param_.max_num_words_per_doc(), param_.max_num_non_zeros_per_word() * param_.word_hashing_param().word_gram_length(), dict_.size() * param_.word_hashing_param().word_gram_length());
        }
    protected:
        /* !
        * \brief For NLettergram, ParseText will return letter grams generated from *UNI_GRAM*
        * This function is only used in building dictionary.
        * Note that when building dictionary, we use UNI_GRAM to make it fair for the first and last word of each document
        */
        virtual vector<string> ParseText(string& str) const override;
        /* !
        * \brief split the token into letter grams
        * \param token to be parsed
        * \return letter grams in a vector
        */
        vector<string> TokenToCharNGram(string& token) const;
    public:
        void Init(const TextHashingParameter& param) override;
        void Extract(const vector<string> &data, SparseBlob<Dtype>& out_blob) const override;
    };
}

#endif