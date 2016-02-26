#pragma once
#include <GDAL/gdal_priv.h>
#include <string>
using namespace std;
class HazePerfection
{
public:
	HazePerfection();
	~HazePerfection();
	void setHotImagePath(const string hotPath);               //设置HOT影像的路径
	void setMaskImagePath(const string maskPath);
	void setMarkImagePath(const string markPath);
	void setTotalChangeImagePath(const string tcPath);
	void setMaxChangeImagePath(const string mcPath);
	void setErodeImagePath(const string erodePath);     //设置形态学侵蚀后的图像路径
	void setResultImagePath(const string resultPath);
	void setErodeTimes(const int times);//设置侵蚀次数
	void setHighThreshold(const float mc, const float tc);//设置maxChange和totalChange的阈值进行二值化
	void setBinaryImagePath(const string binaryPath);

	void morphfill();//填洼算法       1使用setHotImagePath、setInvertA、setMarkB、setMaskImagePath、setMarkImagePath
	void removePeak();  //消峰   2使用setTotalChangeImagePath、setMaxChangeImagePath、setErodeImagePath、setErodeTimes、setHighThreshold、setResultImagePath、setBinaryImagePath


private:
	void readHotImage();                                                    //读取HOT影像
	void invertImage(GDALDataset *pDataset, GDALDataset *dstDataset, float a);                //对图像求反 
	void morphologicalReconstruction(GDALDataset *asMarkDataset, GDALDataset *asMaskDataset);//形态学重建
	void morphologicalErode(GDALDataset *morphDataset);//形态学侵蚀
	void detectionPeak();//检测高亮度区域而非云区域，binaryDataset在此创建

	void setInvertA(const float A);//设置求反时的inverta的值
	void setMarkB(const float B);//设置创建Mark影像时b的值

	int nXSize;
	int nYSize;
	int nBandCount;
	float inverta;//取反时的a得值
	float markB;  //mark影像的值
	int erodeN;  //侵蚀的次数
	float thresholdMC;
	float thresholdTC;
	double sGeoTrans[6];
	void createMask();//创建maskDataset
	void createMark(float b);//创建markDataset
	void createChangeImage();  //创建totalChangeDataset；maxChangeDataset；erodeDataset
	void createResultImage();//创建resultDataset
	float minA(float a[],int n);
	float maxA(float a[], int n);
	void difference(GDALDataset *aDataset, GDALDataset *bDataset);

	string m_hotfilename;
	GDALDataset *hotDataset;   //输入的HOT图像
	string m_maskfilename;
	string m_markfilename;
	GDALDataset *maskDataset;    //填洼算法时的mask影像
	GDALDataset *markDataset;    //填洼算法时的mark影像

	string m_totalchangefilename;
	string m_maxchangefilename;
	string m_erodefileName;
	GDALDataset *totalChangeDataset;    //消峰算法时TC影像
	GDALDataset *maxChangeDataset;     //消峰算法时MC影像
	GDALDataset *erodeDataset;               //形态学侵蚀后的图像

	string m_binaryfilename;
	GDALDataset *binaryDataset;              //满足TC和MC影像的二值化图像，即图像上高亮度而非云区域的影像

	string m_resultfilename;
	GDALDataset *resultDataset;
};

