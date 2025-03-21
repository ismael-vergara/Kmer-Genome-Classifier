#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=LEARN
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Kmer.o \
	${OBJECTDIR}/src/KmerCounter.o \
	${OBJECTDIR}/src/KmerFreq.o \
	${OBJECTDIR}/src/Profile.o \
	${OBJECTDIR}/src/metamain.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-D LEARN -Wall -pedantic
CXXFLAGS=-D LEARN -Wall -pedantic

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CND_CONF}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CND_CONF}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CND_CONF} ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/src/Kmer.o: src/Kmer.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Kmer.o src/Kmer.cpp

${OBJECTDIR}/src/KmerCounter.o: src/KmerCounter.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/KmerCounter.o src/KmerCounter.cpp

${OBJECTDIR}/src/KmerFreq.o: src/KmerFreq.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/KmerFreq.o src/KmerFreq.cpp

${OBJECTDIR}/src/Profile.o: src/Profile.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Profile.o src/Profile.cpp

${OBJECTDIR}/src/metamain.o: src/metamain.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -Iinclude -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/metamain.o src/metamain.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
