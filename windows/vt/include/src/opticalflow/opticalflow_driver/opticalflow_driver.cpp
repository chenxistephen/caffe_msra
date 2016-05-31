#include "..\stdafx.h"
#include "..\cvof2.h"
#include "vtcore.h"

using namespace vt;

int main(int argc, char* argv[])
{
    VT_HR_BEGIN()

    VT_HR_EXIT((argc != 3) ? E_FAIL : S_OK);

    //Read arguments
    WCHAR inputImageName[MAX_PATH] = {0};
    size_t numCharConverted;
    errno_t error = ::mbstowcs_s(&numCharConverted, inputImageName, argv[1], MAX_PATH);
    VT_HR_EXIT((error != 0 || numCharConverted <= 0) ? E_FAIL : S_OK);
    
    // Loads the same image and computes test flow 
    CRGBAImg inputImage1;
    CRGBAImg inputImage2;
    VT_HR_EXIT((VtLoadImage(inputImageName, inputImage1)));
    VT_HR_EXIT((VtLoadImage(inputImageName, inputImage2)));
    
    // Sanity Check Image Sizes Match
    int iWidth1 = inputImage1.Width();
    int iWidth2 = inputImage2.Width();
    int iHeight1 = inputImage1.Height();
    int iHeight2 = inputImage2.Height();
    VT_HR_EXIT((iWidth1 != iWidth2 || iHeight1 != iHeight2) ? E_FAIL : S_OK);

    IOpticalFlow* cHSOpticalFlow;
    
    //Create an optical flow object
    VT_HR_EXIT((vt::CreateOpticalFlow2(cHSOpticalFlow, iWidth1, iHeight1)));
    //Uncomment to create and test the HornSchunck 1-channel flow version
    //VT_HR_EXIT((vt::CreateHornSchunck1Channel(cHSOpticalFlow, iWidth1, iHeight1)));

    // Adds images and computes flow
    VT_HR_EXIT((cHSOpticalFlow->AddImage(inputImage1)));
    VT_HR_EXIT((cHSOpticalFlow->AddImage(inputImage2)));
    VT_HR_EXIT((cHSOpticalFlow->ComputeFlow()));

    //Retrieves flow result and verifies values
    CFloatImg* flowXImage;
    VT_HR_EXIT(cHSOpticalFlow->GetFlowX(flowXImage));
    int flowXWidth = flowXImage->Width();
    int flowXHeight = flowXImage->Height();
    VT_HR_EXIT((flowXWidth != iWidth1 || flowXHeight != iHeight1) ? E_FAIL : S_OK);
    for (int iY = 0; iY < flowXHeight; iY++)
    {
        float* flowXRow = flowXImage->Ptr(iY);
        float* flowXRowEnd = flowXRow + flowXWidth;
        while (flowXRow < flowXRowEnd)
        {
            VT_HR_EXIT((*flowXRow++ == 0) ? S_OK : E_FAIL);
        }
    }

    CFloatImg* flowYImage;
    VT_HR_EXIT(cHSOpticalFlow->GetFlowY(flowYImage));
    int flowYWidth = flowYImage->Width();
    int flowYHeight = flowYImage->Height();
    VT_HR_EXIT((flowYWidth != iWidth1 || flowYHeight != iHeight1) ? E_FAIL : S_OK);
    for (int iY = 0; iY < flowYHeight; iY++)
    {
        float* flowYRow = flowYImage->Ptr(iY);
        float* flowYRowEnd = flowYRow + flowYWidth;
        while (flowYRow < flowYRowEnd)
        {
            VT_HR_EXIT((*flowYRow++ == 0) ? S_OK : E_FAIL);
        }
    }

    VT_HR_EXIT_LABEL()
    if (hr != S_OK)
    {
        ::wprintf(L"Test case failed.\n");
    }
    else
    {
        ::wprintf(L"Test case succeeded.\n");
    }
    return 0;
}