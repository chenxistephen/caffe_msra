#pragma once
#include <memory>

#ifdef IU_SCOPE_PROCESSOR_EXPORTS
#  define IU_WRAPPER_EXT __declspec(dllexport)
#else
#  define IU_WRAPPER_EXT __declspec(dllimport)
#endif 



typedef unsigned char BYTE;
typedef unsigned int UInt32;


namespace IUNativeExecutorWrapper
{
	
	class IU_WRAPPER_EXT IUReadOnlyFeatureWrapper
	{
	public:
		IUReadOnlyFeatureWrapper();
		~IUReadOnlyFeatureWrapper();
	
		class Impl;
		Impl* m_pImpl; 
		
		// these methods provide wrappers for the IIUReadOnlyFeature.
		void Set(const BYTE* pData, UInt32 bufLen, UInt32 width, UInt32 height, UInt32 nChannels);
		void Set(const BYTE* pData, UInt32 bufLen, UInt32 width, UInt32 height, UInt32 nChannels, UInt32 stride);
		const BYTE* GetBuffer() ;
		UInt32 GetWidth() ;
		UInt32 GetHeight() ;
		UInt32 GetStride() ;
		UInt32 GetNumOfChannels() ;
		UInt32 GetLengthInBytes() ;
		UInt32 GetElementCount() ;
		bool IsEmpty() ;
	};



	class IU_WRAPPER_EXT IUNativeExecutor
	{
	public:
		IUNativeExecutor();
		~IUNativeExecutor();      
		bool Initialize(const char* pipelineFilePath, const char **inputFeatureNames, const int inputCount, const char** outputFeatureNames, const int outputCount, const char * dataFolder = ".");

		// The buffers pointed to by outputFeatures[] must be consumed or used by the caller before a subsequent call to Process().
		void Process( IUReadOnlyFeatureWrapper** inputFeatures, int inputCount, IUReadOnlyFeatureWrapper** outputFeatures, int outputCount ) const;	

	private: 
		class Impl;
		Impl* m_pImpl;

		IUNativeExecutor(const IUNativeExecutor &);
		IUNativeExecutor& operator=(const IUNativeExecutor&);

	};
}
