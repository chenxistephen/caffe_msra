// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include <Windows.h>
#include "lib_extract_features.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

using namespace std;

int main(int argc, char** argv)
{
    string modelFilename = "\\CSATRELEVANCE08\\houdong\\src\\MainEmpty\\scratch\\linjuny\\Caffe\\test\\bvlc_reference_caffenet.caffemodel";
    string prototxt = "\\CSATRELEVANCE08\\houdong\\src\\MainEmpty\\scratch\\linjuny\\Caffe\\test\\imagenet_val.prototxt";
    string layerName = "fc7";

    void* Caffe_FE = CaffeCreatFeatureExtractor();
    CaffeLoadModel(Caffe_FE, modelFilename.c_str(), prototxt.c_str());
    

    int dim_features = CaffeGetFeatureDimension(Caffe_FE, layerName.c_str());
    float* outFeature = new float[dim_features];

    ifstream inFile;
    inFile.open("\\CSATRELEVANCE08\\houdong\\src\\MainEmpty\\scratch\\linjuny\\Caffe\\test\\fries1.bgr24");
    stringstream strStream;
    strStream << inFile.rdbuf();
    string str = strStream.str();
    vector<int> strVec;

    boost::char_separator<char> sep(", ");
    boost::tokenizer< boost::char_separator<char> > tokens(str, sep);
    BOOST_FOREACH(const string& t, tokens) {
        strVec.push_back(stoi(t));
    }

    unsigned char* pImageBgrBuffer = new unsigned char[strVec.size()];
    for (int i = 0; i < strVec.size(); i ++){
        pImageBgrBuffer[i] = (unsigned char)strVec[i];
    }
    CaffeExtractFeature(Caffe_FE, pImageBgrBuffer, 240, 186, 720, layerName.c_str(), outFeature);

    std::ofstream myfile("\\CSATRELEVANCE08\\houdong\\src\\MainEmpty\\scratch\\linjuny\\Caffe\\test\\example2.txt");
    for (int d = 0; d < dim_features; ++d) {
        myfile << outFeature[d];
        if (d < dim_features - 1)
            myfile << ",";
    }
    myfile.close();

    return 0;
}