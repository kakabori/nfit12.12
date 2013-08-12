#include <dataset.h>
#include "globalVariables.h"

#include "tvImg.h"

/* Sigma modifications made 11/29/06, but doesn't make any difference
   Search on 'sqrt'  */

int Instruction::qztotnum(){
  /*
    return the total number of qzslices
  */
  int sum=0;
  for(size_t i=0; i<qzstart.size(); i++){
    sum += qzend[i]-qzstart[i]+1;
  }
 printf("qztotsum:%d  \n", sum);
 //printf("qztotsum:%d  backsig:%g \n", sum, backgroundSigmaSquare);
  return sum;
}

int Instruction::qzrange(int &s, int &e){
  /*

  */
  s = qzstart[0];
  e = qzend[0];
  for(size_t i=0; i<qzstart.size(); i++)
    if(qzstart[i] < s)    s = qzstart[i];
  for(size_t i=0; i<qzend.size(); i++)
    if(qzend[i] > e)      e = qzend[i];
  return e-s+1;
}



void Data::exportTiff(char *fname){
  /*
    obsolete, will delete.
  */
  TIFF *tifout;
  uint32 subfiletype=0, imgwidth, imglength;
  uint16 resolutionunit=2;
  uint16 bitspersample=32;
  uint16 compression=1;
  uint16 photometric=1;
  uint16 samplesperpixel=1;
  float xresolution=292.571F,yresolution=292.571F;

  float *buf=NULL;
  unsigned int row,i;
  int s,e;
  imgwidth = ins->qrnum();
  imglength  = ins->qzrange(s,e);
  //printf("%ld %ld\n",imgwidth,imglength);

  if( !(tifout=TIFFOpen(fname,"w")) ) {
    printf("error open tifout\n");
    return;
  }
  TIFFSetField(tifout, TIFFTAG_IMAGEWIDTH, imgwidth);
  TIFFSetField(tifout, TIFFTAG_IMAGELENGTH, imglength);
  TIFFSetField(tifout, TIFFTAG_BITSPERSAMPLE, bitspersample);
  TIFFSetField(tifout, TIFFTAG_DOCUMENTNAME, "can I put something here");
  TIFFSetField(tifout, TIFFTAG_SUBFILETYPE, subfiletype);
  TIFFSetField(tifout, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(tifout, TIFFTAG_PHOTOMETRIC, photometric);
  TIFFSetField(tifout, TIFFTAG_SAMPLESPERPIXEL, samplesperpixel);
  TIFFSetField(tifout, TIFFTAG_ROWSPERSTRIP, imglength);
  TIFFSetField(tifout, TIFFTAG_MINSAMPLEVALUE, 0);
  TIFFSetField(tifout, TIFFTAG_MAXSAMPLEVALUE, 65535);
  TIFFSetField(tifout, TIFFTAG_XRESOLUTION, xresolution);
  TIFFSetField(tifout, TIFFTAG_YRESOLUTION, yresolution);
  TIFFSetField(tifout, TIFFTAG_RESOLUTIONUNIT, resolutionunit);
  TIFFSetField(tifout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  assert(sizeof(float) == 4);
  buf = (float *)malloc( sizeof(float)*imgwidth);
  for (row = 0; row < imglength; row++) {
    for(i=0;i<imgwidth;i++){
      buf[i] = calvalue(i, s+row);
    }
    TIFFWriteScanline(tifout, buf, row, 0);
  }

  TIFFClose(tifout);
  free(buf);
}

void Data::print(FILE *fp){
  /*
    for debugging only. helpper function
    write the dataset to a file
  */
  vector<DataPoint>::iterator it;
  vector<DataPoint>::iterator itend;
  for(size_t i=0; i<qs.size(); i++){
    itend = qs[i].sdata.end();
    for(it = qs[i].sdata.begin(); it != itend; it++){
      fprintf(fp, "%g %g %g\n", it->qx, it->inte, it->cal);
    }
    fprintf(fp,"\n");
  }
}

void Data::writeimg(){
  /*
    update the associated image 'FIT' in %img tog FIT command
	in UserManual
  */
  int width, height, s, e;
  width = ins->qrnum();
  height = ins->qzrange(s,e);

  for(int h=0; h<height; h++){
    for(int w=0; w<width; w++){
      fitimg->setval( w, height-1-h, calvalue(w, s+h) );
    }
  }
}

void Data::writefrm(char *fname, Para *p){
  /*
    output the scaling factors to a file, such as frm.dat
    fname : output filename
    p     : the parameter set used
  */
  FILE *fp = fopen(fname, "w");
  for(size_t i=0; i<qs.size(); i++){
    QzSlice *tqs = &(qs[i]);
    fprintf(fp, "%g %g %g %g %g %g\n", tqs->qz, tqs->scale, tqs->bias,
	    p->setup.getqz(tqs->qz), tqs->sigma_scale, tqs->sigma_bias);
  }
  fclose(fp);
}

void Data::print(char* fn){
  /*
    for debugging only.
    write the dataset to a file
  */
  FILE *fp = fopen(fn, "w");
  print(fp);
  fclose(fp);
}

void Data::info(){
  /*  I/O possibility
    write to standard output the information using
	$dataset info <dataset.name>
	Edit whatever you want to see - very little now.  */

  cout<<"qr:\t"<<ins->qrstart<<'\t'<<ins->qrend<<"\nqz:\t"<<endl;
  for(size_t i=0; i<ins->qzstart.size(); i++){
    cout<<'\t'<<ins->qzstart[i]<<' '<<ins->qzend[i]<<endl;
  }
}

double Data::calvalue(int _qr, int _qz){
  /*
    helper function called by writeimg above.
    Returns calculated value given the pixel position
    note that (_qr, _qz) is the (px,pz) in the thesis

	*/
  int qznum = qs.size();
  for(int k=0;  k<qznum;  k++)
    if( qs[k].qz == _qz )
      return qs[k].sdata[_qr].cal;
  return 0;
}

double Data::datvalue(int _qr, int _qz){
  /* helper -
    return experimental I_e value given the pixel position
    note that (_qr, _qz) is the (px,pz) in the thesis
  */
  int qznum = qs.size();
  for(int k=0;  k<qznum;  k++)
    if( qs[k].qz == _qz )
      return qs[k].sdata[_qr].inte;
  return 0;
}


void Data::readin(YFImgBase *imgptr, Instruction *_ins){
  /*
    initialize the dataset from image pointed by
    'imgptr' with the instructions from '_ins'
  */

  /* bad c++ practice here, will change */
  if(ins != NULL ) delete ins;
  ins = _ins;

  qs.resize(ins->qztotnum());
  n = qs.size() * ins->qrnum();
  int numblocks = ins->qzstart.size();
  int insqrnum = ins->qrnum();
  int i=0;

  for(int k=0;  k<numblocks;  k++){
    int insqzendk = ins->qzend[k];
    for(int m=ins->qzstart[k]; m <= insqzendk; m++, i++){
      qs[i].qz = m;
      qs[i].sdata.resize(insqrnum);
      for(int j=0; j<insqrnum; j++){
	DataPoint *dp;
	dp = &(qs[i].sdata[j]);
	dp->qx = ins->qrstart+j;
	//cout<<m<<' '<<j<<endl; for debugging only
	//cout<<imgptr->getval(m, ins->qrstart+j)<<endl;
	dp->inte = imgptr->getval(ins->qrstart+j, m);
	double avg=0;
	for(int ii=-1; ii<2; ii++){
	  for(int jj=-1; jj<2; jj++){
	    avg += imgptr->getval(ins->qrstart+j+jj, m+ii);
	  }
	}
  /* Old code below.  Data do not justify increase in sigma   xxxx
     with increasing intensity.*/
     //	if(dp->inte < 10) dp->sigma= sqrt( ins->basecnt );
     // else dp->sigma= sqrt( fabs(dp->inte)*(1.0/dupe) + ins->basecnt );
	dp->sigma= sqrt( fabs(dp->inte)*(aFactor/dupe) + ins->basecnt );
      }
    }
  }
  // printf{ " basecnt %g \n " , ins->basecnt ); //debugger xxx
}


void Data::setimg(TVImg<float> *_fitimg){
  /*
    assign an image buffer to the dataset for storing
    the calculated fitting result.
    pre: _fitimg points to an allocated image object
    post: the  image object is resized and positioned.
  */
  fitimg = _fitimg;
  int s, e;
  fitimg->resize( ins->qrnum(), ins->qzrange(s,e) );
  YFDuple<int> rv;
  rv.y = s;
  rv.x = ins->qrstart;
  fitimg->shiftposition(rv);
}

extern double refine(double theor, double exper);

void Data::filter(TVImg<float> *_filterimg){
  /*
    assign an image buffer to the dataset for storing
    the filtered image.
    pre: _filterimg points to an allocated image object
    post: the  image object is resized and positioned
              and filtered inteneities are assigned.
  */
  filterimg = _filterimg;
  int s, e;
  filterimg->resize( ins->qrnum(), ins->qzrange(s,e) );
  YFDuple<int> rv;
  rv.y = s;
  rv.x = ins->qrstart;
  filterimg->shiftposition(rv);

  int width, height;
  width = ins->qrnum();
  height = ins->qzrange(s,e);

  double tmp;
  for(int h=0; h<height; h++){
    for(int w=0; w<width; w++){
      if ( NKfilter!=0 && refine( calvalue(w,s+h), datvalue(w,s+h)) > NKfilter) {tmp=-32000;}
      else {tmp=datvalue(w,s+h);}
      filterimg->setval( w, height-1-h, tmp );
    }
  }
}
