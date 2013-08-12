#include "nfit.h"
#include "modelcalculator.h"
#include "utable.h"
#include <gsl/gsl_sf_bessel.h>
#include <iostream>
#include <fstream>

using namespace std;

#define PI 3.1415926535897932384626433832795
#define SMALLNUM 0.00000000001

void spsumn(ModelCalculator& mc, double qz)
{
	mc.setSliceParameter(qz);
	mc.QxSlice(0, 0.302);
	
	vector<double> r;
	vector<double> sumn;
	mc.get_spsumnPoints(r, sumn);
	
	vector<double>::iterator iter, iter2;
	iter2 = sumn.begin();
	int count = 0;
	for (iter = r.begin(); iter != r.end(); iter++) {
		cout << *iter << " " << exp(*iter) - SMALLNUM << " " << *iter2 << endl;
		iter2++;
		count++;
	}
	cout << count << endl;
}

void rotated(ModelCalculator& mc, double qz)
{
	vector<double> qxvec;
	vector<double> rotatedSFvec;
	mc.getRotatedSF(0.008, 0.3, qz, qxvec, rotatedSFvec);
	
	ofstream myfile;
	myfile.open("rotatedSF_qz0.45_new.dat");
	for (size_t i = 0; i < rotatedSFvec.size(); i++) {
		myfile << qxvec[i] << " " << rotatedSFvec[i] << endl;
		//cout << qxvec[i] << " " << rotatedSFvec[i] << endl;
	}
	myfile.close();

	vector<double> qx, rotatedSF;
	mc.get_spRotatedPoints(qx, rotatedSF);
	vector<double>::iterator iter, iter2;
	iter2 = rotatedSF.begin();
	int count = 0;
	for (iter = qx.begin(); iter != qx.end(); iter++) {
		cout << *iter << " " << exp(*iter) - SMALLNUM << " " << exp(*iter2) << endl;
		iter2++;
		count++;
	}
	cout << count << endl;

}

void StrFct(ModelCalculator& mc, double qz)
{
	mc.setSliceParameter(qz);
	mc.QxSlice(0.00, 0.3);
	vector<double> qvec, sfvec;
	mc.get_spStrFctPoints(qvec, sfvec);
	cout << qvec.size() << endl;
	for (vector<double>::size_type i = 0; i < qvec.size(); i++) {
		cout << exp(qvec[i]) - SMALLNUM << " " << exp(sfvec[i]) << endl;
	}
}

void CCD(ModelCalculator& mc, double qz)
{
	mc.setSliceParameter(qz);
	mc.QxSlice(0.00, 0.3);
	
	for (double qx = 0.00; qx < 0.3;){
		cout << qx << " " << mc.getCCDStrFct(qx) << endl;
		qx += 0.01;
	}
}

int main()
{
	//double t; int n;
	Utable utab;
	utab.writeUtableFile("utab_nfit12.15.dat", 600, 2000);
	utab.readUtableFile("utab_nfit12.15.dat");
	//cout << "Enter n and t: ";
	//while (cin >> n >> t) {
	//	cout << utab.interp_utable(n, t) << " "
	//	     << utab.gsl_interp_utable(n, t) << endl;
	//}
	
	ModelCalculator mc;
	mc.read_in_utable("utab_nfit12.15.dat");
	  
	mc.setModelParameter(5.6212e-13, "Kc");
	mc.setModelParameter(4.82716e13, "B");
	mc.setModelParameter(60.519, "D");
	mc.setModelParameter(37, "T");
	mc.setModelParameter(10000, "Lr");
	mc.setModelParameter(6.35298, "Mz");
	mc.setModelParameter(0.0134, "edisp");
	mc.setModelParameter(0, "mosaic");
	mc.setModelParameter(1.175, "wavelength");
	mc.setModelParameter(0.07113, "pixelSize");
	mc.setModelParameter(2.3, "bFWHM");
	mc.setModelParameter(359.3, "s");

	//double qz = 0.534023;
	double qz = 0.5;
	//CCD(mc, qz);

	StrFct(mc, qz);
	//rotated(mc, qz);
	//spsumn(mc, qz);	
	//spHr(mc);
	//bessel();
	//create2DMap(mc);
	//cout << getLowerLimit(0.1, 0.415, 1.175) << endl;
	//cout << getUpperLimit(0.1, 0.415, 1.175) << endl;

	return 0;
}

