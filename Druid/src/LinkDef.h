#ifndef EVENTDISPLAY_LINKDEF_H_
#define EVENTDISPLAY_LINKDEF_H_

#ifdef __CINT__

/** $Id: EventDisplay_LinkDef.h,v 1.2 2007-03-13 22:50:24 sanchez Exp $
 * 
 * LinkDef file needed to generate ROOT dictionary for EventDisplay.
 * To do so:
 * rootcint -f EventDisplayDict.cc -c ../include/EventDisplay.h EventDisplay_LinkDef.h
 * 
 * @ author   Allister Levi Sanchez (LLR)
 * @ date     2007.02.26
 */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//#pragma link C++ namespace XED;
#pragma link C++ nestedclasses;
//#pragma link C++ class XED::GeomGearILD00;
//#pragma link C++ class XED::IGeometry;
#pragma link C++ class EventNavigator;


#endif

#endif /*EVENTDISPLAY_LINKDEF_H_*/

