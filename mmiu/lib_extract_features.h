#ifdef __CAFFE_DLL_EXPORT
#define CAFFE_DLL_EXPORT __declspec(dllexport)
#else
#define CAFFE_DLL_EXPORT __declspec(dllimport)
#endif

#include<iostream>


CAFFE_DLL_EXPORT class IFeatureExtractor
{
public:
	virtual bool LoadModel(const char* modelFilename, const char* prototxt) = 0;
    virtual void ExtractCaffeFeature(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature) = 0;
    virtual void ExtractCaffeImageFeature(unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature) = 0;
    virtual void ExtractCaffeTextFeature(unsigned char* pTextBuffer, const int length, const char* layerName, float* outFeature) = 0;
    virtual void ExtractCaffeVectorFeature(unsigned char* pVectorBuffer, const int length, const char* layerName, float* outFeature) = 0;
    virtual size_t GetFeatureDimension(const char* layerName) = 0;
    virtual void Destroy() = 0;
};


extern "C"
{
	CAFFE_DLL_EXPORT void* CaffeCreatFeatureExtractor();
	CAFFE_DLL_EXPORT bool CaffeLoadModel(void* pFeatureExtractor, const char* modelFilename, const char* prototxt);
    CAFFE_DLL_EXPORT void CaffeExtractFeature(void* pFeatureExtractor, unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature);
    CAFFE_DLL_EXPORT void CaffeExtractImageFeature(void* pFeatureExtractor, unsigned char* pImageBgrBuffer, const int width, const int height, const int stride, const char* layerName, float* outFeature);
    CAFFE_DLL_EXPORT void CaffeExtractTextFeature(void* pFeatureExtractor, unsigned char* pTextBuffer, const int length, const char* layerName, float* outFeature);
    CAFFE_DLL_EXPORT void CaffeExtractVectorFeature(void* pFeatureExtractor, unsigned char* pVectorBuffer, const int length, const char* layerName, float* outFeature);
	CAFFE_DLL_EXPORT size_t CaffeGetFeatureDimension(void* pFeatureExtractor, const char* layerName);
	CAFFE_DLL_EXPORT void CaffeDestroy(void* pFeatureExtractor);
}