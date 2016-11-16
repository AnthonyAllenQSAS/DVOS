
# Mac OSX
DOT_SO=.dylib
OSX_RELEASE=OSX
LD=g++ -stdlib=libstdc++  -O3  -dynamiclib -fPIC
CPL=g++ -stdlib=libstdc++ -fPIC -I. -Werror 

# Linux
#DOT_SO=.so
#LD=g++ -std=c++11'  -O3  -dynamiclib 
#CPL=g++ -std=c++11' -I. -Werror

# Windows
#DOT_SO=.dll
#LD=/mingw64/bin/g++  -shared  -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc -mthreads -Wl,-subsystem,windows -Wl,--stack,800000000
#CPL=/mingw64/bin/g++ -O3 -fno-keep-inline-dllexport -Wno-deprecated -mthreads -DWINDOWS 


TARGET_LIB=../lib

.cc.o:
	${CPL} -c $*.cc


DVOS_OBJS=\
DvObj.o \
DvTime.o \
DvOperators.o \
DvUtilities.o \
DvMatrixOps.o \
DvAlgorithms.o \
DvUnit.o \
DvJoin.o \
DvComponent.o \
DvString.o \
DvXML.o


all: clean new

clean:
	'rm' -f  *.o
	'rm' -f  *.a
	'rm' -f  *.la
	'rm' -f  *.so
	'rm' -f *.dylib
	'rm' -f *.dll
	'rm' -f *~
	'rm' -f *%
	'rm' -f ../lib/libdvos.so
	'rm' -f ../lib/libdvos.dylib
	'rm' -f ../lib/libdvos.dll

new:  libdvos.so

libdvos.so: install $(DVOS_OBJS)
	${LD} -o libdvos$(DOT_SO) $(DVOS_OBJS) 
	cp  libdvos$(DOT_SO) ${TARGET_LIB}
ifdef OSX_RELEASE
	install_name_tool -id "${TARGET_LIB}/libdvos${DOT_SO}" ${TARGET_LIB}/libdvos${DOT_SO}
endif


install:
	cp DvObj.h ../include
	cp Dvar.h ../include
	cp DvTime.h ../include
	cp DvString.h ../include
	cp DvRecord.h ../include
	cp DvUnit.h ../include
	cp Xrefs.h ../include
 