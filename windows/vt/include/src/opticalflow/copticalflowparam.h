//+---------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation 2010.  All rights reserved.
//
//  Description:
//     Parameter Class for Optical Flow Algorithms
//
//  History:
//      2011/8/12-wkienzle/sbaker
//          Created
//
//----------------------------------------------------------------------------

#pragma once
#include "windows.h"

/// <summary> An parameter object for IOpticalFlow optical flow algorithms. All 
/// parameters are represented as strings. Parameters have an optional 
/// description field. </summary>
class COFParam
{
public:     
	 /// <summary> Constructor - sets all strings to be empty </summary>
     COFParam()
      {
            m_name[0] = '\0';
            m_value[0] = '\0';
            m_description[0] = '\0';
      }

	  /// <summary> Sets the name string of the parameter </summary>
      HRESULT SetName(const char* name)
      {
            return Set(m_name, name);
      }

 	  /// <summary> Gets the name string of the parameter </summary>
     const char* GetName() const
      {
            return m_name;
      }

	  /// <summary> Sets the value string of the parameter </summary>
      HRESULT SetValue(const char* value)
      {
            return Set(m_value, value);
      }

	  /// <summary> Gets the value string of the parameter </summary>
      const char* GetValue() const
      {
            return m_value;
      }

	  /// <summary> Sets the description string of the parameter </summary>
      HRESULT SetDescription(const char* description)
      {
            return Set(m_description, description);
      }

 	  /// <summary> Gets the description string of the parameter </summary>
      const char* GetDescription() const
      {
            return m_description;
      }

	  /// <summary> The maximum length of the name string </summary>
      static const size_t MAX_NAME_LENGTH = 256;
	  /// <summary> The maximum length of the value string </summary>
      static const size_t MAX_VALUE_LENGTH = 256;
	  /// <summary> The maximum length of the description string </summary>
      static const size_t MAX_DESCRIPTION_LENGTH = 256;
      
private:
	  /// <summary> The name string </summary>
      char m_name[MAX_NAME_LENGTH + 1];
 	  /// <summary> The value string </summary>
      char m_value[MAX_VALUE_LENGTH + 1];
	  /// <summary> The description string </summary>
      char m_description[MAX_DESCRIPTION_LENGTH + 1];

 	  /// <summary> Function to set strings </summary>
      template<size_t ArraySize>
      static HRESULT Set(char (&dst)[ArraySize], const char* src)
      {
            if (strlen(src) > ArraySize - 1)
            {
				return E_INVALIDARG;
			}
            if (strcpy_s(dst, src) != 0) 
            {
				return E_FAIL;
			}
            return S_OK;
      }
};

