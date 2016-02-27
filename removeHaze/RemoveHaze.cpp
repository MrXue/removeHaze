#include "RemoveHaze.h"
#include <vector>
#include <algorithm>

HazePerfection::HazePerfection() :nXSize(0), nYSize(0), nBandCount(1), inverta(0.0), markB(0.0), 
erodeN(1), thresholdMC(0.0), thresholdTC(0.0), filterTemplateN(3)
{
}


HazePerfection::~HazePerfection()
{
	if (hotDataset != NULL)
		GDALClose((GDALDatasetH)hotDataset);
	if (maskDataset != NULL)
		GDALClose((GDALDatasetH)maskDataset);
	if (markDataset!=NULL)
		GDALClose((GDALDatasetH)markDataset);
	if (totalChangeDataset != NULL)
		GDALClose((GDALDatasetH)totalChangeDataset);
	if (maxChangeDataset != NULL)
		GDALClose((GDALDatasetH)maxChangeDataset);
	if (erodeDataset != NULL)
		GDALClose((GDALDatasetH)erodeDataset);
	if (binaryDataset != NULL)
		GDALClose((GDALDatasetH)binaryDataset);
	if (resultDataset != NULL)
		GDALClose((GDALDatasetH)resultDataset);
}

void HazePerfection::setHotImagePath(const string hotPath)
{
	m_hotfilename = hotPath;
}

void HazePerfection::readHotImage()
{
	if (m_hotfilename.empty())
		return;
	const char *hotfile = m_hotfilename.c_str();
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	hotDataset = (GDALDataset *)GDALOpen(hotfile, GA_ReadOnly);
	if (hotDataset == nullptr)
		return;
	nXSize = hotDataset->GetRasterXSize();
	nYSize = hotDataset->GetRasterYSize();
	hotDataset->GetGeoTransform(sGeoTrans);
	double minAndmax[2];
	hotDataset->GetRasterBand(1)->ComputeRasterMinMax(false, minAndmax);
	setInvertA(minAndmax[1] + 5.0);
	setMarkB(0.0);
}

void HazePerfection::invertImage(GDALDataset *pDataset, GDALDataset *dstDataset, float a)
{
	GDALRasterBand *band = pDataset->GetRasterBand(1);
	GDALDataType dataType = band->GetRasterDataType();
	float *pixelData = new float[nXSize*nYSize];
	band->RasterIO(GF_Read, 0, 0, nXSize, nYSize, pixelData, nXSize, nYSize, GDT_Float32, 0, 0);
	
	GDALRasterBand *dstBand = dstDataset->GetRasterBand(1);
	for (int j = 0; j < nYSize; j++)
	{
		for (int i = 0; i < nXSize; i++)
		{
			pixelData[j*nXSize + i] = a - pixelData[j*nXSize + i];
		}
	}
	dstBand->RasterIO(GF_Write, 0, 0, nXSize, nYSize, pixelData, nXSize, nYSize, GDT_Float32, 0, 0);
	delete[]pixelData;
}

void HazePerfection::createMask()
{
	const char *pszFormat = "GTiff";
	const char *maskfile = m_maskfilename.c_str();
	GDALDriver *poDriver = (GDALDriver*)GDALGetDriverByName(pszFormat);
	maskDataset = poDriver->Create(maskfile, nXSize, nYSize, 1, GDT_Float32, NULL);
	maskDataset->SetGeoTransform(sGeoTrans);
	maskDataset->SetProjection(hotDataset->GetProjectionRef());
}

void HazePerfection::createMark(float b)
{
	const char *pszFormat = "GTiff";
	const char *markfile = m_markfilename.c_str();
	GDALDriver *poDriver = (GDALDriver*)GDALGetDriverByName(pszFormat);
	markDataset = poDriver->Create(markfile, nXSize, nYSize, 1, GDT_Float32, NULL);
	markDataset->SetGeoTransform(sGeoTrans);
	markDataset->SetProjection(hotDataset->GetProjectionRef());
	float *pixelData = new float[nXSize*nYSize];
	GDALRasterBand *maskBand = maskDataset->GetRasterBand(1);
	maskBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize, pixelData, nXSize, nYSize, GDT_Float32, 0, 0);
	for (int i = 1; i < nYSize-1; i++)
	{
		for (int j = 1; j < nXSize - 1; j++)
			pixelData[i*nXSize + j] = b;
	}
	GDALRasterBand *markBand = markDataset->GetRasterBand(1);
	markBand->RasterIO(GF_Write, 0, 0, nXSize, nYSize, pixelData, nXSize, nYSize, GDT_Float32, 0, 0);
	delete[]pixelData;
}

void HazePerfection::morphologicalReconstruction(GDALDataset *asMarkDataset, GDALDataset *asMaskDataset)
{
	int nmarkXSize = asMarkDataset->GetRasterXSize();
	int nmarkYSize = asMarkDataset->GetRasterYSize();
	int nmaskXSize = asMaskDataset->GetRasterXSize();
	int nmaskYSize = asMaskDataset->GetRasterYSize();
	//Mark图像和Mask图像维数不一样
	if (nmaskXSize != nmarkXSize || nmarkYSize != nmaskYSize)
		return;
	float *markPixelData = new float[nmarkXSize*nmarkYSize];
	float *maskPixelData = new float[nmaskXSize*nmaskYSize];
	asMarkDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nmarkXSize, nmarkYSize, markPixelData, nmarkXSize, nmarkYSize, GDT_Float32, 0, 0);
	asMaskDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nmaskXSize, nmaskYSize, maskPixelData, nmaskXSize, nmaskYSize, GDT_Float32, 0, 0);
	float temp[3] = { 0, 0, 0 };
	for (int i = 1; i < nmarkYSize; i++)
	{
		for (int j = 1; j < nmarkXSize; j++)
		{
			//取大
			temp[0] = markPixelData[i*nmarkXSize + j - 1];
			temp[1] = markPixelData[(i - 1)*nmarkXSize + j];
			temp[2] = markPixelData[i*nmarkXSize + j];
			markPixelData[i*nmarkXSize + j] = maxA(temp, 3);
			/*markPixelData[i*nmarkXSize + j - 1] >= markPixelData[(i - 1)*nmarkXSize + j] ? markPixelData[i*nmarkXSize + j] = markPixelData[i*nmarkXSize + j - 1] : markPixelData[i*nmarkXSize + j] = markPixelData[(i - 1)*nmarkXSize + j];
			if (markPixelData[i*nmarkXSize + j]<markPixelData[(i - 1)*nmarkXSize + j - 1])
			markPixelData[i*nmarkXSize + j] = markPixelData[(i - 1)*nmarkXSize + j - 1];*/
			//取小
			if (markPixelData[i*nmarkXSize + j] > maskPixelData[i*nmaskXSize + j])
				markPixelData[i*nmarkXSize + j] = maskPixelData[i*nmaskXSize + j];
		}
		for (int j = nmarkXSize - 2; j>=0; j--)
		{
			//取大
			temp[0] = markPixelData[i*nmarkXSize + j + 1];
			temp[1] = markPixelData[(i - 1)*nmarkXSize + j];
			temp[2] = markPixelData[i*nmarkXSize + j];
			markPixelData[i*nmarkXSize + j] = maxA(temp,3);
			/*markPixelData[i*nmarkXSize + j + 1] >= markPixelData[(i - 1)*nmarkXSize + j] ? markPixelData[i*nmarkXSize + j] = markPixelData[i*nmarkXSize + j + 1] : markPixelData[i*nmarkXSize + j] = markPixelData[(i - 1)*nmarkXSize + j];
			if (markPixelData[i*nmarkXSize + j]<markPixelData[(i - 1)*nmarkXSize + j + 1])
			markPixelData[i*nmarkXSize + j] = markPixelData[(i - 1)*nmarkXSize + j + 1];*/
			//取小
			if (markPixelData[i*nmarkXSize + j] > maskPixelData[i*nmaskXSize + j])
				markPixelData[i*nmarkXSize + j] = maskPixelData[i*nmaskXSize + j];
		}
	}

	for (int i = nmarkYSize - 2; i >= 0; i--)
	{
		for (int j = 1; j < nmarkXSize; j++)
		{
			temp[0] = markPixelData[i*nmarkXSize + j - 1];
			temp[1] = markPixelData[(i + 1)*nmarkXSize + j];
			temp[2] = markPixelData[i*nmarkXSize + j];
			markPixelData[i*nmarkXSize + j] = maxA(temp,3);
			/*markPixelData[i*nmarkXSize + j - 1] >= markPixelData[(i + 1)*nmarkXSize + j] ? markPixelData[i*nmarkXSize + j] = markPixelData[i*nmarkXSize + j - 1] : markPixelData[i*nmarkXSize + j] = markPixelData[(i + 1)*nmarkXSize + j];
			if (markPixelData[i*nmarkXSize + j]<markPixelData[(i + 1)*nmarkXSize + j - 1])
				markPixelData[i*nmarkXSize + j] = markPixelData[(i + 1)*nmarkXSize + j - 1];*/
			if (markPixelData[i*nmarkXSize + j]>maskPixelData[i*nmaskXSize + j])
				markPixelData[i*nmarkXSize + j] = maskPixelData[i*nmaskXSize + j];
		}
		for (int j = nmarkXSize - 2; j >= 0; j--)
		{
			temp[0] = markPixelData[i*nmarkXSize + j + 1];
			temp[1] = markPixelData[(i + 1)*nmarkXSize + j];
			temp[2] = markPixelData[i*nmarkXSize + j];
			markPixelData[i*nmarkXSize + j] = maxA(temp,3);
			/*markPixelData[i*nmarkXSize + j + 1] >= markPixelData[(i + 1)*nmarkXSize + j] ? markPixelData[i*nmarkXSize + j] = markPixelData[i*nmarkXSize + j + 1] : markPixelData[i*nmarkXSize + j] = markPixelData[(i + 1)*nmarkXSize + j];
			if (markPixelData[i*nmarkXSize + j]<markPixelData[(i + 1)*nmarkXSize + j + 1])
			markPixelData[i*nmarkXSize + j] = markPixelData[(i + 1)*nmarkXSize + j + 1];*/
			if (markPixelData[i*nmarkXSize + j] > maskPixelData[i*nmaskXSize + j])
				markPixelData[i*nmarkXSize + j] = maskPixelData[i*nmaskXSize + j];
		}
	}
	asMarkDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, nmarkXSize, nmarkYSize, markPixelData, nmarkXSize, nmarkYSize, GDT_Float32, 0, 0);
	delete[]markPixelData;
	delete[]maskPixelData;
}

void HazePerfection::morphfill()
{
	readHotImage();
	createMask();
	invertImage(hotDataset, maskDataset, inverta);
	createMark(markB);
	morphologicalReconstruction(markDataset, maskDataset);
	invertImage(markDataset, maskDataset, inverta);
}

void HazePerfection::setInvertA(const float A)
{
	inverta = A;
}

void HazePerfection::setMarkB(const float B)
{
	markB = B;
}

void HazePerfection::morphologicalErode(GDALDataset *morphDataset)
{
	createChangeImage();
	GDALRasterBand *morphBand = morphDataset->GetRasterBand(1);
	float *morphPixelsData = new float[nXSize*nYSize];
	float *maxPixelsData = new float[nXSize*nYSize];
	float *erodePixelsData = new float[nXSize*nYSize];
	float *tempErodeData = new float[nXSize*nYSize];
	float tempChange = 0.0;
	
	float a[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	morphBand->RasterIO(GF_Read, 0, 0, nXSize, nYSize, morphPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	for (int k = 0; k < erodeN; k++)
	{
		if (k == 0){
			for (int i = 1; i < nYSize - 1; i++)
			{
				for (int j = 1; j < nXSize - 1; j++)
				{
					a[0] = morphPixelsData[(i - 1)*nXSize + j - 1];
					a[1] = morphPixelsData[(i - 1)*nXSize + j];
					a[2] = morphPixelsData[(i - 1)*nXSize + j + 1];
					a[3] = morphPixelsData[i*nXSize + j - 1];
					a[4] = morphPixelsData[i*nXSize + j];
					a[5] = morphPixelsData[i*nXSize + j + 1];
					a[6] = morphPixelsData[(i + 1)*nXSize + j - 1];
					a[7] = morphPixelsData[(i + 1)*nXSize + j];
					a[8] = morphPixelsData[(i + 1)*nXSize + j + 1];
					tempErodeData[i*nXSize + j] = minA(a, 9);
					erodePixelsData[i*nXSize + j] = tempErodeData[i*nXSize + j];
					maxPixelsData[i*nXSize + j] = morphPixelsData[i*nXSize + j] - erodePixelsData[i*nXSize + j];
				}
			}
			for (int i = 0; i < nXSize; i++)
			{
				tempErodeData[i] = morphPixelsData[i];
				tempErodeData[nXSize*(nYSize - 1) + i] = morphPixelsData[nXSize*(nYSize - 1) + i];
				erodePixelsData[i] = morphPixelsData[i];
				erodePixelsData[nXSize*(nYSize - 1) + i] = morphPixelsData[nXSize*(nYSize - 1) + i];
				maxPixelsData[i] = 0.0;
				maxPixelsData[nXSize*(nYSize - 1) + i] = 0.0;
			}
			/*tempErodeData[0] = morphPixelsData[0];
			erodePixelsData[0] = morphPixelsData[0];
			maxPixelsData[0] = 0.0;
			tempErodeData[nXSize*nYSize - 1] = morphPixelsData[nXSize*nYSize];
			erodePixelsData[nXSize*nYSize - 1] = morphPixelsData[nXSize*nYSize];
			maxPixelsData[nXSize*nYSize - 1] = 0.0;*/
			for (int j = 1; j < nYSize; j++)
			{
				tempErodeData[j*nXSize] = morphPixelsData[j*nXSize];
				tempErodeData[j*nXSize - 1] = morphPixelsData[j*nXSize - 1];
				erodePixelsData[j*nXSize] = morphPixelsData[j*nXSize];
				erodePixelsData[j*nXSize - 1] = morphPixelsData[j*nXSize - 1];
				maxPixelsData[j*nXSize] = 0.0;
				maxPixelsData[j*nXSize - 1] = 0.0;
			}
		}
		else
		{
			for (int i = 1; i < nYSize - 1;i++)
			{
				for (int j = 1; j < nXSize - 1; j++)
				{
					a[0] = erodePixelsData[(i - 1)*nXSize + j - 1];
					a[1] = erodePixelsData[(i - 1)*nXSize + j];
					a[2] = erodePixelsData[(i - 1)*nXSize + j + 1];
					a[3] = erodePixelsData[i*nXSize + j - 1];
					a[4] = erodePixelsData[i*nXSize + j];
					a[5] = erodePixelsData[i*nXSize + j + 1];
					a[6] = erodePixelsData[(i + 1)*nXSize + j - 1];
					a[7] = erodePixelsData[(i + 1)*nXSize + j];
					a[8] = erodePixelsData[(i + 1)*nXSize + j + 1];
					erodePixelsData[i*nXSize + j] = minA(a, 9);
					tempChange = tempErodeData[i*nXSize + j] - erodePixelsData[i*nXSize + j];
					tempErodeData[i*nXSize + j] = erodePixelsData[i*nXSize + j];
					if (maxPixelsData[i*nXSize + j] < tempChange)
						maxPixelsData[i*nXSize + j] = tempChange;
				}
			}
		}
	}
	GDALRasterBand *erodeBand = erodeDataset->GetRasterBand(1);
	GDALRasterBand *maxBand = maxChangeDataset->GetRasterBand(1);
	/*for (int i = 0; i < nYSize; i++)
	{
		for (int j = 0; j < nXSize; j++)
		{
			tempErodeData[i*nXSize + j] = morphPixelsData[i*nXSize + j] - erodePixelsData[i*nXSize + j];
		}
	}*/
	erodeBand->RasterIO(GF_Write, 0, 0, nXSize, nYSize, erodePixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	maxBand->RasterIO(GF_Write, 0, 0, nXSize, nYSize, maxPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	//totalChangeDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, nXSize, nYSize, tempErodeData, nXSize, nYSize, GDT_Float32, 0, 0);
	delete[]erodePixelsData;
	delete[]maxPixelsData;
	delete[]morphPixelsData;
	delete[]tempErodeData;
}

void HazePerfection::setErodeTimes(const int times)
{
	erodeN = times;
}

float HazePerfection::minA(float a[], int n)
{
	float temp = a[0];
	for (int i = 1; i < n; i++)
	{
		if (a[i] < temp)
			temp = a[i];
	}
	return temp;
}

void HazePerfection::createChangeImage()
{
	const char *pszFormat = "GTiff";
	const char *totalFile = m_totalchangefilename.c_str();
	GDALDriver *poDrive = (GDALDriver*)GDALGetDriverByName(pszFormat);
	totalChangeDataset = poDrive->Create(totalFile, nXSize, nYSize, 1, GDT_Float32, NULL);
	totalChangeDataset->SetGeoTransform(sGeoTrans);
	totalChangeDataset->SetProjection(hotDataset->GetProjectionRef());

	const char *maxfile = m_maxchangefilename.c_str();
	GDALDriver *pomaxDrive = (GDALDriver*)GDALGetDriverByName(pszFormat);
	maxChangeDataset = pomaxDrive->Create(maxfile, nXSize, nYSize,1, GDT_Float32,NULL);
	maxChangeDataset->SetGeoTransform(sGeoTrans);
	maxChangeDataset->SetProjection(hotDataset->GetProjectionRef());

	const char *erodefile = m_erodefileName.c_str();
	GDALDriver *erodeDrive = (GDALDriver*)GDALGetDriverByName(pszFormat);
	erodeDataset = erodeDrive->Create(erodefile, nXSize, nYSize, 1, GDT_Float32, NULL);
	erodeDataset->SetGeoTransform(sGeoTrans);
	erodeDataset->SetProjection(hotDataset->GetProjectionRef());
}


void HazePerfection::setMaskImagePath(const string maskPath)
{
	m_maskfilename = maskPath;
}

void HazePerfection::setMarkImagePath(const string markPath)
{
	m_markfilename = markPath;
}

void HazePerfection::setTotalChangeImagePath(const string tcPath)
{
	m_totalchangefilename = tcPath;
}

void HazePerfection::setMaxChangeImagePath(const string mcPath)
{
	m_maxchangefilename = mcPath;
}

void HazePerfection::setErodeImagePath(const string erodePath)
{
	m_erodefileName = erodePath;
}

void HazePerfection::detectionPeak()
{
	morphologicalErode(maskDataset);
	morphologicalReconstruction(erodeDataset, maskDataset);
	difference(maskDataset, erodeDataset);
	const char *pszFormat = "GTiff";
	const char *binaryfile = m_binaryfilename.c_str();
	GDALDriver *poDrive = (GDALDriver*)GDALGetDriverByName(pszFormat);
	binaryDataset = poDrive->Create(binaryfile, nXSize, nYSize, 1, GDT_Float32, NULL);
	binaryDataset->SetGeoTransform(sGeoTrans);
	binaryDataset->SetProjection(hotDataset->GetProjectionRef());

	float *binaryPixelsData = new float[nXSize*nYSize];
	float *tempPixelsData = new float[nXSize*nYSize];
	totalChangeDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nXSize, nYSize, binaryPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	maxChangeDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nXSize, nYSize, tempPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	for (int i = 0; i < nYSize; i++)
	{
		for (int j = 0; j < nXSize; j++)
		{
			if (binaryPixelsData[i*nXSize + j]>=thresholdTC && tempPixelsData[i*nXSize + j]>=thresholdMC)
				binaryPixelsData[i*nXSize + j] = 1.0;
			else
			{
				binaryPixelsData[i*nXSize + j] = 0.0;
			}
		}
	}
	binaryDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, nXSize, nYSize, binaryPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	delete[]binaryPixelsData;
	delete[]tempPixelsData;
}

void HazePerfection::setHighThreshold(const float mc, const float tc)
{
	thresholdMC=mc;
	thresholdTC = tc;
}

//erodeDataset;binaryDataset
void HazePerfection::removePeak()
{
	detectionPeak();
	createResultImage();
	//float a[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	float *erodePixelsData = new float[nXSize*nYSize];
	float *binaryPixelsData = new float[nXSize*nYSize];
	maskDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nXSize, nYSize, erodePixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	binaryDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nXSize, nYSize, binaryPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	int t = (int)(filterTemplateN / 2);
	int temp = t;
	int n = t;
	for (int i = t; i < nYSize - t; i++)
	{
		for (int j = t; j < nXSize - t; j++)
		{
			if (binaryPixelsData[i*nXSize + j] == 1.0)
			{
				std::vector<float>filter;
				while (t != 0)
				{
					while (temp != 0)
					{
						filter.push_back(erodePixelsData[(i - t)*nXSize + j - temp]);
						filter.push_back(erodePixelsData[(i - t)*nXSize + j + temp]);
						filter.push_back(erodePixelsData[(i + t)*nXSize + j - temp]);
						filter.push_back(erodePixelsData[(i + t)*nXSize + j + temp]);
						temp--;
					}
					temp = n;
					t--;
				}
				t = n;
				while (t != 0)
				{
					filter.push_back(erodePixelsData[i*nXSize + j + t]);
					filter.push_back(erodePixelsData[i*nXSize + j - t]);
					t--;
				}
				t = n;
				filter.push_back(erodePixelsData[i*nXSize + j]);
				std::sort(filter.begin(), filter.end());
				binaryPixelsData[i*nXSize + j] = *filter.begin();
			}
			else
			{
				binaryPixelsData[i*nXSize + j] = erodePixelsData[i*nXSize + j];
			}
		}
	}

	////3*3最小值滤波
	//for (int i = 1; i < nYSize-1; i++)
	//{
	//	for (int j = 1; j < nXSize-1; j++)
	//	{
	//		if (binaryPixelsData[i*nXSize + j] == 1.0)
	//		{
	//			a[0] = erodePixelsData[(i - 1)*nXSize + j - 1];
	//			a[1] = erodePixelsData[(i - 1)*nXSize + j ];
	//			a[2] = erodePixelsData[(i - 1)*nXSize + j + 1];
	//			a[3] = erodePixelsData[i*nXSize + j - 1];
	//			a[4] = erodePixelsData[i*nXSize + j ];
	//			a[5] = erodePixelsData[i*nXSize + j + 1];
	//			a[6] = erodePixelsData[(i + 1)*nXSize + j - 1];
	//			a[7] = erodePixelsData[(i + 1)*nXSize + j ];
	//			a[8] = erodePixelsData[(i + 1)*nXSize + j + 1];
	//			binaryPixelsData[i*nXSize + j] = minA(a, 9);
	//		}
	//		else
	//		{
	//			binaryPixelsData[i*nXSize + j] = erodePixelsData[i*nXSize + j];
	//		}
	//	}
	//}
	resultDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, nXSize, nYSize, binaryPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	delete[]erodePixelsData;
	delete[]binaryPixelsData;
}

void HazePerfection::createResultImage()
{
	const char *pszFormat = "GTiff";
	const char *resultfile = m_resultfilename.c_str();
	GDALDriver *poDriver = (GDALDriver *)GDALGetDriverByName(pszFormat);
	resultDataset = poDriver->Create(resultfile, nXSize, nYSize, 1, GDT_Float32, NULL);
	resultDataset->SetGeoTransform(sGeoTrans);
	resultDataset->SetProjection(hotDataset->GetProjectionRef());
}

void HazePerfection::setResultImagePath(const string resultPath)
{
	m_resultfilename = resultPath;
}

void HazePerfection::setBinaryImagePath(const string binaryPath)
{
	m_binaryfilename = binaryPath;
}

float HazePerfection::maxA(float a[], int n)
{
	float temp = a[0];
	for (int i = 1; i < n; i++)
	{
		if (temp < a[i])
			temp = a[i];
	}
	return temp;
}

void HazePerfection::difference(GDALDataset *aDataset, GDALDataset *bDataset)
{
	float *aPixelsData = new float[nXSize*nYSize];
	float *bPixelsData = new float[nXSize*nYSize];
	aDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nXSize, nYSize, aPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	bDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, nXSize, nYSize, bPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	for (int i = 0; i < nYSize; i++)
	{
		for (int j = 0; j < nXSize; j++)
			aPixelsData[i*nXSize + j] = aPixelsData[i*nXSize + j] - bPixelsData[i*nXSize + j];
	}
	totalChangeDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, nXSize, nYSize, aPixelsData, nXSize, nYSize, GDT_Float32, 0, 0);
	delete[]aPixelsData;
	delete[]bPixelsData;
}

void HazePerfection::setFilterTemplate(int n /*= 3*/)
{
	int temp = n % 2;
	if (temp == 0)
		filterTemplateN = 3;
	else
	{
		filterTemplateN = n;
	}
}

GDALDataset * HazePerfection::getResultDataset()
{
	return resultDataset;
}
