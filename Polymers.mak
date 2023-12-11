# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=Polymers - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to Polymers - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "Polymers - Win32 Release" && "$(CFG)" !=\
 "Polymers - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "Polymers.mak" CFG="Polymers - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Polymers - Win32 Release" (based on\
 "Win32 (x86) Console Application")
!MESSAGE "Polymers - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "Polymers - Win32 Debug"
F90=fl32.exe
RSC=rc.exe

!IF  "$(CFG)" == "Polymers - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\Polymers.exe"

CLEAN : 
	-@erase ".\Release\Polymers.exe"
	-@erase ".\Release\Input.obj"
	-@erase ".\Release\FI_POL.obj"
	-@erase ".\Release\DICOTOMY.obj"
	-@erase ".\Release\MK_DO_MAIN.obj"
	-@erase ".\Release\Numerical_Fmixing_calculation.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Polymers.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/Polymers.pdb" /machine:I386 /out:"$(OUTDIR)/Polymers.exe" 
LINK32_OBJS= \
	"$(INTDIR)/Input.obj" \
	"$(INTDIR)/FI_POL.obj" \
	"$(INTDIR)/DICOTOMY.obj" \
	"$(INTDIR)/MK_DO_MAIN.obj" \
	"$(INTDIR)/Numerical_Fmixing_calculation.obj" \
	"..\..\LIB\MATHS.LIB" \
	"..\..\LIB\MATHD.LIB"

"$(OUTDIR)\Polymers.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "Polymers - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\Polymers.exe"

CLEAN : 
	-@erase ".\Debug\Polymers.exe"
	-@erase ".\Debug\FI_POL.obj"
	-@erase ".\Debug\Input.obj"
	-@erase ".\Debug\DICOTOMY.obj"
	-@erase ".\Debug\MK_DO_MAIN.obj"
	-@erase ".\Debug\Numerical_Fmixing_calculation.obj"
	-@erase ".\Debug\Polymers.ilk"
	-@erase ".\Debug\Polymers.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Zi /I "Debug/" /c /nologo
F90_PROJ=/Zi /I "Debug/" /c /nologo /Fo"Debug/" /Fd"Debug/Polymers.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/Polymers.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/Polymers.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/Polymers.exe" 
LINK32_OBJS= \
	"$(INTDIR)/FI_POL.obj" \
	"$(INTDIR)/Input.obj" \
	"$(INTDIR)/DICOTOMY.obj" \
	"$(INTDIR)/MK_DO_MAIN.obj" \
	"$(INTDIR)/Numerical_Fmixing_calculation.obj" \
	"..\..\LIB\MATHS.LIB" \
	"..\..\LIB\MATHD.LIB"

"$(OUTDIR)\Polymers.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "Polymers - Win32 Release"
# Name "Polymers - Win32 Debug"

!IF  "$(CFG)" == "Polymers - Win32 Release"

!ELSEIF  "$(CFG)" == "Polymers - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\MK_DO_MAIN.f90

"$(INTDIR)\MK_DO_MAIN.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\Input.f90

"$(INTDIR)\Input.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\FI_POL.f90

"$(INTDIR)\FI_POL.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\DICOTOMY.f90

"$(INTDIR)\DICOTOMY.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=\MSDEV\LIB\MATHS.LIB

!IF  "$(CFG)" == "Polymers - Win32 Release"

!ELSEIF  "$(CFG)" == "Polymers - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=\MSDEV\LIB\MATHD.LIB

!IF  "$(CFG)" == "Polymers - Win32 Release"

!ELSEIF  "$(CFG)" == "Polymers - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\Numerical_Fmixing_calculation.f90

"$(INTDIR)\Numerical_Fmixing_calculation.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
