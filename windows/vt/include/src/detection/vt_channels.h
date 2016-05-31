#pragma once

#include "vtcore.h"
#include "vt_planarimage.h"

namespace vt_internal
{
	template <typename StaticParams>
	class NormalizationEngine;
	template <typename StaticParams>
	class GradientEngine;
	template <typename StaticParams>
	class BinningEngine;
}

namespace vt
{
	struct CChannelsDefaultStaticParams
	{
		struct InputNormParams
		{
			static const unsigned int normFactorLUTPreShift = 0;
			static const unsigned int numNormFactorLUTInputBits = 8;
			static const unsigned int normFactorFixShift = 13;
			static const unsigned int log2NumNormFactorLUTBins = 8;
			typedef unsigned short NormFactorLUTBinType;
		};

		struct GradientParams
		{
			static const unsigned int numMagLUTInputBits = 8;
			static const unsigned int log2NumMagLUTBins = 7;
			static const unsigned int numMagBits = 8;
			typedef unsigned char MagLUTBinType;

			static const unsigned int numAngLUTInputBits = 8;
			static const unsigned int log2NumAngLUTBins = 7;
			static const unsigned int numAngBits = 8;
			typedef unsigned char AngLUTBinType;
		};

		struct GradientNormParams
		{
			static const unsigned int normFactorLUTPreShift = 0;
			static const unsigned int numNormFactorLUTInputBits = 8;
			static const unsigned int normFactorFixShift = 13;
			static const unsigned int log2NumNormFactorLUTBins = 8;
			typedef unsigned short NormFactorLUTBinType;
		};

		struct BinningParams
		{
			static const unsigned int log2OrientBinSize = 2;
			static const unsigned int numOrientBinLUTInputBits = 8;
			static const unsigned int log2NumOrientBinLUTBins = 8;
			typedef int OrientLUTBinType;
			static const unsigned int orientBinWeightFixShift = 4;
		};
	};

	struct CChannelsParams
	{		
		CChannelsParams()
			: inputIsTransposed(false),
			numInputChannels(1),
			shrinkFactor(4),
			numInputBlurIterations(1),
			numInputNormBlurIterations(7),
			inputNormRegularization(25),
			inputNormScaling(0.75f),
			numGradNormBlurIterations(2),
			gradNormRegularization(15),
			gradNormScaling(0.45f),
			numOrientBins(6)
		{
		}

		// input orientation (for Matlab compatibility)
		bool inputIsTransposed;

		// number of input channels
		unsigned int numInputChannels;

		// input to channel downsampling
		unsigned int shrinkFactor;

		// input blur
		unsigned int numInputBlurIterations;

		// input normalization
		unsigned int numInputNormBlurIterations;
		unsigned int inputNormRegularization;
		float inputNormScaling;

		// gradient normalization
		unsigned int numGradNormBlurIterations;
		unsigned int gradNormRegularization;
		float gradNormScaling;

		// binning
		unsigned int numOrientBins;
	};

	enum ChannelType
	{
		ChannelTypeUnknown,
		ChannelTypeInput,
		ChannelTypeGradient
	};

	template <typename StaticParams>
	class CChannelsEngine
	{
	public:
		CChannelsEngine();
		virtual ~CChannelsEngine();

	private:
		CChannelsEngine(const CChannelsEngine&);
		CChannelsEngine& operator=(const CChannelsEngine&);

	public:
		HRESULT Initialize(unsigned int maxInputWidth, 
			unsigned int maxInputHeight, const CChannelsParams& params);

		bool IsInitialized() const;
		HRESULT GetParams(CChannelsParams& p) const;

		HRESULT GetNumChannels(unsigned int& numChannels) const;
		HRESULT GetChannelType(ChannelType& type, unsigned int index) const;

		HRESULT GetChannelSize(CSize& channelSize, unsigned int inputWidth, 
			unsigned int inputHeight) const;		
		HRESULT GetInputSize(CSize& inputSize, unsigned int channelsWidth, 
			unsigned int channelsHeight) const;

		HRESULT ComputeChannels(CPlanarByteImg& channels, 
			const CByteImg& input) const;

	private:
		static bool ValidateParams(const CChannelsParams& params);

		CChannelsParams params;	

		unsigned int maxInputWidth;
		unsigned int maxInputHeight;

		bool isInitialized;		

	private:
		static unsigned int GetNumChannels(const CChannelsParams& params);

		bool GetInputBuffer(CPlanarByteImg& inByteBuf, 
			unsigned int inputWidth, unsigned int inputHeight) const;

		mutable CPlanarByteImg inputByteBuffers;
		mutable CCompositeImg<LumaType<int>> channelRowIntBuffers;

	private:
		typedef typename StaticParams::InputNormParams StaticInputNormParams;
		typedef typename StaticParams::GradientParams StaticGradientParams;
		typedef typename StaticParams::GradientNormParams 
			StaticGradientNormParams;
		typedef typename StaticParams::BinningParams StaticBinningParams;

		vt_internal::NormalizationEngine<StaticInputNormParams>* 
			inputNormEngine;
		vt_internal::GradientEngine<StaticGradientParams>* 
			gradEngine;
		vt_internal::NormalizationEngine<
			StaticGradientNormParams>* gradNormEngine;
		vt_internal::BinningEngine<StaticBinningParams>* 
			binEngine;
	};
}