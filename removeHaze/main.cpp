#include <iostream>
#include "RemoveHaze.h"
#include <string>
using namespace std;

int main()
{
	string hotImagePath;
	cout << "输入HOT影像的路径：" << endl;
	cin >> hotImagePath;
	string maskImagePath;
	cout << "创建mask影像的路径：" << endl;
	cin >> maskImagePath;
	string markImagePath;
	cout << "创建mark影像的路径：" << endl;
	cin >> markImagePath;
	string tcImagePath;
	cout << "创建TotalChange影像的路径：" << endl;
	cin >> tcImagePath;
	string mcImagePath;
	cout << "创建MaxChange影像的路径：" << endl;
	cin >> mcImagePath;
	string erodeIamgePath;
	cout << "创建侵蚀后影像的路径：" << endl;
	cin >> erodeIamgePath;
	string binaryImagePath;
	cout << "创建二值化后影像的路径：" << endl;
	cin >> binaryImagePath;
	string resultImagePath;
	cout << "创建结果影像的路径：" << endl;
	cin >> resultImagePath;
	HazePerfection *testHazePer = new HazePerfection;
	testHazePer->setHotImagePath(hotImagePath);
	testHazePer->setMaskImagePath(maskImagePath);
	testHazePer->setMarkImagePath(markImagePath);
	testHazePer->setTotalChangeImagePath(tcImagePath);
	testHazePer->setMaxChangeImagePath(mcImagePath);
	testHazePer->setErodeImagePath(erodeIamgePath);
	testHazePer->setBinaryImagePath(binaryImagePath);
	testHazePer->setResultImagePath(resultImagePath);

	//TODO:以下值待确定
	testHazePer->setErodeTimes(5);
	testHazePer->setHighThreshold(10.0, 10.0);

	testHazePer->morphfill();
	testHazePer->removePeak();
	return 1;
}