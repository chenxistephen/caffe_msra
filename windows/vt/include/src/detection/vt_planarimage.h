#pragma once

#include "vtcore.h"

namespace vt 
{
#ifdef MAKE_DOXYGEN_WORK
	void foo();
#endif

	template <typename T>
	class CPlanarTypedImg
	{
	public:
		typedef CCompositeImg<LumaType<T>,T> CPlaneType;

	public:
		/// <summary> Constructor </summary> 
		CPlanarTypedImg()
			: m_numBands(0), m_height(0)
		{
		}

		/// <summary> Destructor </summary> 
		virtual ~CPlanarTypedImg()
		{
		}

	private:
		// copying is verboten	
		CPlanarTypedImg(const CPlanarTypedImg&);
		CPlanarTypedImg &operator=(const CPlanarTypedImg&);

	public:
		/// <summary> Initialize image, allocate pixel memory. </summary>
		/// <returns> 
		///		- S_OK on success
		///		- E_INVALIDARG if requested properties are invalid
		///		- E_OUTOFMEMORY if memory allocation failed 
		/// </returns>
		HRESULT Create(int iW, int iH, int iBands, 
			AlignMode eAlign = DEFAULT_ALIGN)
		{ 
			VT_HR_BEGIN();

			m_height = iH;
			m_numBands = iBands;

			VT_HR_EXIT(m_storage.Create(iW, iH * iBands, eAlign));

			VT_HR_EXIT(CreateBandViews());

			VT_HR_END();
		}

		/// <summary> Initialize image, wrap existing pixel buffer. </summary>
		/// <returns> 
		///		- S_OK on success
		///		- E_INVALIDARG if requested properties are invalid or if 
		/// pbBuffer is NULL.
		/// </returns>
		/// <DL><DT> Remarks: </DT></DL>
		///		- No memory is allocated
		///		- Note that the pixel buffer is not const, since it may be 
		///	modified later via this CPlanarTypedImg.
		///		- The wrapped memory resource is managed by the caller. 
		HRESULT Create(T *pbBuffer, int iW, int iH, int iBands, 
			int iStrideBytes)
		{
			VT_HR_BEGIN();

			m_height = iH;
			m_numBands = iBands;

			VT_HR_EXIT(m_storage.Create(reinterpret_cast<Byte*>(pbBuffer), 
				iW, iH * iBands, iStrideBytes));

			VT_HR_EXIT(CreateBandViews());

			VT_HR_END();
		}

		/// <summary> Share pixels. </summary>
		///	<returns> 
		///		- true, if the image has been shared successfully
		/// </returns>
		HRESULT Share(CPlanarTypedImg<T>& other) const
		{
			return Share(other, 0, Bands());
		}

		/// <summary> Share pixels (with band selection). </summary>
		///	<returns> 
		///		- true, if the image has been shared successfully
		/// </returns>
		HRESULT Share(CPlanarTypedImg<T>& other, int firstBand, int numBands) 
			const
		{
			VT_HR_BEGIN();

			if (firstBand < 0 || firstBand >= Bands())
				return E_INVALIDARG;
			if (numBands <= 0 || firstBand + numBands - 1 >= Bands())
				return E_INVALIDARG;

			other.m_height = m_height;
			other.m_numBands = numBands;
			
			CRect shareRect(0, Height() * firstBand, Width(), 
				Height() * (firstBand + numBands));

			VT_HR_EXIT(m_storage.Share(other.m_storage, &shareRect));

			VT_HR_EXIT(other.CreateBandViews());

			VT_HR_END();
		}
		
		/// <summary> Return if image is valid. </summary>
		///	<returns> 
		///		- true, if the image has been created successfully
		/// </returns>
		bool IsValid() const
		{
			return m_storage.IsValid();
		}

		/// <summary> Return image width in pixels. </summary>
		int Width() const
		{
			return m_storage.Width();
		}

		/// <summary> Return number of bytes per row. </summary>
		int StrideBytes() const
		{
			return m_storage.StrideBytes();
		}

		/// <summary> Return image height in pixels. </summary>
		int Height() const
		{
			return m_height;
		}

		/// <summary> Return number of bands. </summary>
		int Bands() const
		{
			return m_numBands;
		}

		/// <summary> Return const reference to band iBand. 
		/// WARNING #1: if the image is not initialized, or iBand is out of
		/// range, calling Band() will lead to undefined behavior. 
		/// WARNING #2: the returned reference points to an internal storage 
		/// that is not aware this reference's existence. Thus, the returned 
		/// reference becomes invalid if *this changes or is destroyed. Also, 
		/// altering the storage of the returned object (e.g. by calling 
		/// Create()) will corrupt *this CTypedPlanarImage and result in 
		/// undefined behavior. 
		/// </summary>
		const CPlaneType& Band(int iBand) const
		{
			VT_ASSERT(iBand >= 0 && iBand < m_numBands);

			return m_bands[iBand];
		}

		/// <summary> Return reference to band iBand (see const overload for 
		/// detailed description). </summary>
		CPlaneType& Band(int iBand)
		{
			return const_cast<CPlaneType&>(
				static_cast<const CPlanarTypedImg&>(*this).Band(iBand));
		}

	private:
		HRESULT CreateBandViews()
		{
			VT_HR_BEGIN();

			VT_HR_EXIT(m_bands.resize(m_numBands));

			for (int i = 0; i < m_numBands; ++i)
			{
				VT_HR_EXIT(m_bands[i].Create(
					const_cast<Byte*>(m_storage.BytePtr(i * Height())), 
					Width(), Height(), StrideBytes()));
			}

			VT_HR_END();
		}

	private:
		int m_height;
		int m_numBands;

		CPlaneType m_storage;
		vector<CPlaneType> m_bands;
	};

	typedef CPlanarTypedImg<Byte> CPlanarByteImg;
	typedef CPlanarTypedImg<UInt16> CPlanarShortImg;
	typedef CPlanarTypedImg<float> CPlanarFloatImg;
	typedef CPlanarTypedImg<HALF_FLOAT> CPlanarHalfFloatImg;

	template <typename OT, typename IT>
	HRESULT VtConvertPlanarToPackedImage(CTypedImg<OT>& dst, 
		const CPlanarTypedImg<IT>& src)
	{
		VT_HR_BEGIN();

		VT_HR_EXIT(src.IsValid() ? S_OK : E_INVALIDSRC);		

		if (dst.Width() != src.Width() || dst.Height() != src.Height() 
			|| dst.Bands() != src.Bands())
		{
			VT_HR_EXIT(dst.Create(src.Width(), src.Height(), src.Bands()));
		}

		if (src.Bands() == 1)
		{
			VT_HR_EXIT(VtConvertImage(dst, src.Band(0)));
			return S_OK;
		}

		for (int y=0; y<src.Height(); ++y)
		{
			for (int x=0; x<src.Width(); ++x)
			{
				for (int b=0; b<src.Bands(); ++b)
				{
					VtConv(dst.Ptr(x, y, b), src.Band(b).Pix(x, y));
				}
			}
		}

		VT_HR_END();
	}

	template <typename OT, typename IT>
	HRESULT VtConvertPackedToPlanarImage(CPlanarTypedImg<OT>& dst, 
		const CTypedImg<IT>& src)
	{
		VT_HR_BEGIN();

		VT_HR_EXIT(src.IsValid() ? S_OK : E_INVALIDSRC);		

		if (dst.Width() != src.Width() || dst.Height() != src.Height() 
			|| dst.Bands() != src.Bands())
		{
			VT_HR_EXIT(dst.Create(src.Width(), src.Height(), src.Bands()));
		}

		if (src.Bands() == 1)
		{
			VT_HR_EXIT(VtConvertImage(dst.Band(0), src));
			return S_OK;
		}

		for (int y=0; y<src.Height(); ++y)
		{
			for (int x=0; x<src.Width(); ++x)
			{
				for (int b=0; b<src.Bands(); ++b)
				{
					VtConv(dst.Band(b).Ptr(x, y), src.Pix(x, y, b));
				}
			}
		}

		VT_HR_END();
	}

	template <typename OT, typename IT>
	HRESULT VtConvertPlanarImage(CPlanarTypedImg<OT>& dst, 
		const CPlanarTypedImg<IT>& src)
	{
		VT_HR_BEGIN();

		VT_HR_EXIT(src.IsValid() ? S_OK : E_INVALIDSRC);		

		if (dst.Width() != src.Width() || dst.Height() != src.Height() 
			|| dst.Bands() != src.Bands())
		{
			VT_HR_EXIT(dst.Create(src.Width(), src.Height(), src.Bands()));
		}

		for (int b = 0; b < src.Bands(); ++b)
			VT_HR_EXIT(VtConvertImage(dst.Band(b), src.Band(b)));
		
		VT_HR_END();
	}
};