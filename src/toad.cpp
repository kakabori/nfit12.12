/*  This is the MAIN program  */

#include <cstdlib>
#include <cstring>
#include <tcl.h>
#include <tk.h>
#include <blt.h>
#include <tiff.h>
#include <tiffio.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include "tvImg.h"
#include "toad.h"
#include "toadmisc.h"
#include "toadcmd.h"

extern "C" int Toad_Init(Tcl_Interp *interp);
extern GlobalData gData;

/* This starts an event loop in Tcl that never returns
   Tcl then takes all commands from the console  */

int main(int argc, char *argv[]){
  printf( " Start main  ");
  gData.interp=Tcl_CreateInterp();
  Tk_MainEx(argc, argv, Toad_Init,gData.interp);
  return(0);
}
