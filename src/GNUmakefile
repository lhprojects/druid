#  $Id: GNUmakefile,v 1.15 2007-03-26 13:22:24 sanchez Exp $
#
#  Makefile for ILC Event Display, a C++ program for viewing an event
#  at the ILC detector / CALICE testbeam by specifying the event number and
#  the GEAR XML file (which describes the detector geometry)
#
#  2007.03.12  Allister Levi Sanchez (LLR)

PROGRAM = Druid
EXE = ../bin/$(PROGRAM)

CC = g++ $(CFLAGS)
LD = $(CC)

ROOTINC = `root-config --cflags`
LCIOINC = -I$(LCIO)/include
#GEARINC = -I$(GEAR)/include

CFLAGS = -g -O -I../include $(ROOTINC) $(LCIOINC) $(DPMINC)

ROOTLIBS = `root-config --glibs` -lGeom -lGeomPainter -lGeomBuilder -lGed
ROOTLIBS    += -lEve -lTreePlayer -lGeom -lRGL -lGed -lEG

LCIOLIBS = -I $(wildcard ${LCIO}/include) -L $(wildcard ${LCIO}/lib*) -llcio -lsio -lz
#LCIOLIBS = -Wl,-rpath,$(LCIO)/lib64 -lsio -lz 
#GEARLIBS = -L$(GEAR)/lib -lgearxml -lgear
LDFLAGS = $(ROOTLIBS) $(LCIOLIBS) $(GEARLIBS)

DICTSRCS = GlobalDict.cc
DICTHDRS = $(patsubst %.cc,%.h,$(DICTSRCS))

SRCS = BuildGeo.cc CaloHits.cc EventNavigator.cc GUI.cc loadevent.cc BuildPFOs.cc ClusterHits.cc Vertex.cc\
	BuildMCParticles.cc Druid.cc TrackerHits.cc TrackAssignedHits.cc MultiView.cc Connect.cc RecoJets.cc\
        geometry.cc trajectory.cc segment3.cc point.cc fitting_root.cc\
        $(DICTSRCS)

OBJS = $(patsubst %.cc,%.o,$(SRCS))

$(info "info: LDFLAGS = $(LDFLAGS)")

all: $(OBJS)
	if [ ! -d ../bin ]; then mkdir ../bin; fi;
	$(LD) -o $(EXE) $(OBJS) $(LDFLAGS) -Wl,-rpath,$(wildcard $(LCIO)/lib*)
	cp GlobalDict_rdict.pcm ../bin/GlobalDict_rdict.pcm

$(DICTSRCS): ../include/*.hh ../include/*.h 
	rootcint -f $@ -c $^ $(LCIOINC) LinkDef.h

.cc.o: $(SRCS) 
	$(CC) -c $(patsubst %.o,%.cc,$@)

clean:
	rm -f $(EXE) $(OBJS) *~ ../include/*~ $(DICTSRCS) $(DICTHDRS)
	rm -rf ../bin
	rm -f GlobalDict_rdict.pcm

