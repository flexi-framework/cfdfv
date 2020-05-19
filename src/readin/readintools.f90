MODULE MOD_ReadInTools
!===================================================================================================================================
! ReadinTools 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE ISO_VARYING_STRING
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC::GETSTR
PUBLIC::CNTSTR
PUBLIC::GETINT
PUBLIC::GETREAL
PUBLIC::GETLOGICAL
PUBLIC::GETINTARRAY
PUBLIC::GETREALARRAY
PUBLIC::IgnoredStrings
!===================================================================================================================================
INTERFACE GETSTR
  MODULE PROCEDURE GETSTR
END INTERFACE

INTERFACE CNTSTR
  MODULE PROCEDURE CNTSTR
END INTERFACE

INTERFACE GETINT
  MODULE PROCEDURE GETINT
END INTERFACE

INTERFACE GETREAL
  MODULE PROCEDURE GETREAL
END INTERFACE

INTERFACE GETLOGICAL
  MODULE PROCEDURE GETLOGICAL
END INTERFACE

INTERFACE GETINTARRAY
  MODULE PROCEDURE GETINTARRAY
END INTERFACE

INTERFACE GETREALARRAY
  MODULE PROCEDURE GETREALARRAY
END INTERFACE

INTERFACE IgnoredStrings
  MODULE PROCEDURE IgnoredStrings
END INTERFACE

INTERFACE FillStrings
  MODULE PROCEDURE FillStrings
END INTERFACE

INTERFACE FindStr
  MODULE PROCEDURE FindStr
END INTERFACE

INTERFACE LowCase
  MODULE PROCEDURE LowCase
END INTERFACE

INTERFACE GetNewString
  MODULE PROCEDURE GetNewString
END INTERFACE

INTERFACE DeleteString
  MODULE PROCEDURE DeleteString
END INTERFACE

TYPE tString
  TYPE(Varying_String)::Str
  TYPE(tString),POINTER::NextStr,PrevStr
END TYPE tString

LOGICAL,PUBLIC::ReadInDone=.FALSE.
TYPE(tString),POINTER::FirstString

CONTAINS

FUNCTION GETSTR(Key,Proposal)
!===================================================================================================================================
! Read string named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255)                   :: GetStr   ! String read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,GetStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,GetStr,DefMsg)
END IF
WRITE(UNIT_StdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM(Key),' | ', TRIM(GetStr),' | ',TRIM(DefMsg),' | '
END FUNCTION GETSTR



FUNCTION CNTSTR(Key,Proposal)
!===================================================================================================================================
! Counts all occurances of string named "key" from inifile and store in "GETSTR". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETSTR" (error if "Proposal" not given).
! Inifile was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
INTEGER,OPTIONAL,INTENT(IN)          :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: CntStr   ! Number of parameters named "Key" in inifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=LEN(Key))              :: TmpKey
TYPE(tString),POINTER                :: Str1
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

CntStr=0
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)

! Search
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) THEN
    CntStr=CntStr+1
  END IF ! (INDEX...
  ! Next string in list
  Str1=>Str1%NextStr
END DO
IF (CntStr.EQ.0) THEN
  IF (PRESENT(Proposal)) THEN
    CntStr=Proposal
  ELSE
    WRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
    STOP
  END IF
END IF
END FUNCTION CNTSTR



FUNCTION GETINT(Key,Proposal)
!===================================================================================================================================
! Read integer named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: GetInt  ! Integer read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetInt
WRITE(UNIT_StdOut,'(a3,a30,a3,i33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetInt,' | ',TRIM(DefMsg),' | '
END FUNCTION GETINT



FUNCTION GETREAL(Key,Proposal)
!===================================================================================================================================
! Read real named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: GetReal  ! Real read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetReal
WRITE(UNIT_StdOut,'(a3,a30,a3,e33.5,a3,a7,a3)')' | ',TRIM(Key),' | ', GetReal,' | ',TRIM(DefMsg),' | '
END FUNCTION GETREAL



FUNCTION GETLOGICAL(Key,Proposal)
!===================================================================================================================================
! Read logical named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key        ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal   ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: GetLogical ! Logical read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr
CHARACTER(LEN=8)                     :: DefMsg
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetLogical
WRITE(UNIT_StdOut,'(a3,a30,a3,l33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetLogical,' | ',TRIM(DefMsg),' | '
END FUNCTION GETLOGICAL



FUNCTION GETINTARRAY(Key,nIntegers,Proposal)
!===================================================================================================================================
! Read array of "nIntegers" integer values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
! list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nIntegers        ! Number of values in array
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal         ! Default values as character string (as in setup file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                   :: GetIntArray(nIntegers)      ! Integer array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)        :: HelpStr
CHARACTER(LEN=8)          :: DefMsg
INTEGER                   :: iInteger
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetIntArray
WRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                               'Integer array of size (',nIntegers,') | ',TRIM(DefMsg),' | '
DO iInteger=0,nIntegers-1
  IF ((iInteger.GT.0) .AND. (MOD(iInteger,8).EQ.0)) THEN
    WRITE(UNIT_stdOut,*)
    WRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
  END IF
  WRITE(UNIT_stdOut,'(i5)',ADVANCE='NO')GetIntArray(iInteger+1)
END DO
WRITE(UNIT_stdOut,*)
END FUNCTION GETINTARRAY



FUNCTION GETREALARRAY(Key,nReals,Proposal)
!===================================================================================================================================
! Read array of "nReals" real values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
! list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nReals           ! Number of values in array
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal         ! Default values as character string (as in setup file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                      :: GetRealArray(nReals)        ! Real array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)        :: HelpStr
CHARACTER(LEN=8)          :: DefMsg
INTEGER                   :: iReal
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings
ReadInDone=.TRUE.

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*)GetRealArray
WRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                               'Real array of size (',nReals,') | ',TRIM(DefMsg),' | '
  WRITE(UNIT_stdOut,'(A)') 
DO iReal=0,nReals-1
  IF ((iReal.GT.0) .AND. (MOD(iReal,8).EQ.0)) THEN
    WRITE(UNIT_stdOut,*)
    WRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
  END IF
  WRITE(UNIT_stdOut,'(f12.6)',ADVANCE='NO')GetRealArray(iReal+1)
END DO
WRITE(UNIT_stdOut,*)
END FUNCTION GETREALARRAY



SUBROUTINE IgnoredStrings()
!===================================================================================================================================
! Prints out remaining strings in list after read-in is complete
!===================================================================================================================================
! MODULES
USE ISO_VARYING_STRING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tString),POINTER                  :: Str1
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)')" THE FOLLOWING INI-FILE PARAMETERS WERE IGNORED:"
Str1=>FirstString
DO WHILE(ASSOCIATED(Str1))
  WRITE(UNIT_stdOut,'(A4,A)')" |- ",TRIM(CHAR(Str1%Str))
  Str1=>Str1%NextStr
END DO
WRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE IgnoredStrings



SUBROUTINE FillStrings(IniFile)
!===================================================================================================================================
! Read ini file and put each line in a string object. All string objects are connected to a list of string objects starting
! with "firstString"
!===================================================================================================================================
! MODULES
USE ISO_VARYING_STRING
USE,INTRINSIC :: ISO_FORTRAN_ENV
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL   :: IniFile                    ! Name of ini file to be read in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tString),POINTER                  :: Str1=>NULL(),Str2=>NULL()
CHARACTER(LEN=255)                     :: HelpStr,Str
CHARACTER(LEN=300)                     :: File
TYPE(Varying_String)                   :: aStr,bStr,Separator
INTEGER                                :: EOF
!===================================================================================================================================
! Check if we have read in ini file already
IF (ReadInDone) RETURN
! Get name of ini file
IF (PRESENT(IniFile)) THEN
  File = TRIM(IniFile)
ELSE
  CALL GETARG(1,File)
  !CALL GET_COMMAND_ARGUMENT(1,File)
END IF
WRITE(UNIT_StdOut,*)'| Reading from file "',TRIM(File),'":'

OPEN(UNIT   = 103,        &
     FILE   = File,       &
     STATUS = 'OLD',      &
     ACTION = 'READ',     &
     ACCESS = 'SEQUENTIAL')
EOF=0

NULLIFY(Str1,Str2)
DO WHILE(EOF.NE.IOSTAT_END)
  IF(.NOT.ASSOCIATED(Str1)) CALL GetNewString(Str1)
    ! Read line from file
    CALL Get(103,aStr,iostat=EOF)
    Str=aStr
!IPWRITE(*,*)'Reading: ',Str,EOF
    IF (EOF.NE.IOSTAT_END) THEN
      ! Remove comments with "!"
      CALL Split(aStr,Str1%Str,"!")
      ! Remove comments with "#"
      CALL Split(Str1%Str,bStr,"#")
      Str1%Str=bStr
      ! Remove "%" sign from old ini files, i.e. mesh% disc% etc.
      CALL Split(Str1%Str,bStr,"%",Separator,Back=.false.)
      ! If we have a newtype ini file, take the other part
      IF(LEN(CHAR(Separator)).EQ.0) Str1%Str=bStr
      ! Remove blanks
      Str1%Str=Replace(Str1%Str," ","",Every=.true.)
      ! Replace brackets
      Str1%Str=Replace(Str1%Str,"(/"," ",Every=.true.)
      Str1%Str=Replace(Str1%Str,"/)"," ",Every=.true.)
      ! Replace commas
      Str1%Str=Replace(Str1%Str,","," ",Every=.true.)
      ! Lower case
      CALL LowCase(CHAR(Str1%Str),HelpStr)
      ! If we have a remainder (no comment only)
      IF(LEN_TRIM(HelpStr).GT.2) THEN
        Str1%Str=Var_Str(HelpStr)
        IF(.NOT.ASSOCIATED(Str2)) THEN
          FirstString=>Str1
        ELSE
          Str2%NextStr=>Str1
          Str1%PrevStr=>Str2
        END IF
        Str2=>Str1
        CALL GetNewString(Str1)
      END IF
    END IF
END DO
CLOSE(103)
END SUBROUTINE FillStrings



SUBROUTINE GetNewString(Str)
!===================================================================================================================================
! Create and initialize new string object.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tString),POINTER :: Str ! New string
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
NULLIFY(Str)
ALLOCATE(Str)
NULLIFY(Str%NextStr,Str%PrevStr)
END SUBROUTINE GetNewString



SUBROUTINE DeleteString(Str)
!===================================================================================================================================
! Remove string "Str" from list of strings witFirstString,h first element "DirstString" and delete string.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    TYPE(tString),POINTER :: Str         ! String to delete
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF (ASSOCIATED(Str%NextStr)) Str%NextStr%PrevStr=>Str%PrevStr
IF (ASSOCIATED(Str,FirstString)) THEN
  FirstString=>Str%NextStr
ELSE
  Str%PrevStr%NextStr=>Str%NextStr
END IF
DEALLOCATE(Str)
NULLIFY(Str)
END SUBROUTINE DeleteString



SUBROUTINE FindStr(Key,Str,DefMsg,Proposal)
!===================================================================================================================================
! Find parameter string containing keyword "Key" in list of strings starting with "FirstString" and return string "Str" without
! keyword. If keyword is not found in list of strings, return default values "Proposal" (error if not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key         ! Search for this keyword in ini file
CHARACTER(LEN=8),INTENT(INOUT)       :: DefMsg      ! Default message = keyword not found, return default parameters (if available)
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal    ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT)         :: Str         ! Parameter string without keyword
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=LEN(Key))              :: TmpKey 
TYPE(tString),POINTER                :: Str1
LOGICAL                              :: Found
!===================================================================================================================================
DefMsg='*CUSTOM'
! Convert to lower case
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)
Found=.FALSE.
Str1=>FirstString
DO WHILE(.NOT.Found)
  IF (.NOT.ASSOCIATED(Str1)) THEN
    IF (.NOT.PRESENT(Proposal)) THEN
      WRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
      STOP
    ELSE ! Return default value
      CALL LowCase(TRIM(Proposal),Str)
      IF (Str(1:1).NE.'@') THEN
        DefMsg='DEFAULT'
      END IF
      RETURN
    END IF ! (.NOT.PRESENT(Proposal))
  END IF ! (.NOT.ASSOCIATED(Str1))

  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) THEN
    Found=.TRUE.
    Str1%Str=replace(Str1%Str,TRIM(TmpKey)//'=',"",Every=.TRUE.)
    Str=TRIM(CHAR(Str1%Str))
    ! Remove string from list
    CALL DeleteString(Str1)
  ELSE
    ! Next string in list
    Str1=>Str1%NextStr
  END IF

END DO
END SUBROUTINE FindStr



SUBROUTINE LowCase(Str1,Str2)
!===================================================================================================================================
! Transform upper case letters in "Str1" into lower case letters, result is "Str2"
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: Str1 ! Input string 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT) :: Str2 ! Output string, lower case letters only
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: iLen,nLen,Upper
CHARACTER(LEN=*),PARAMETER   :: lc='abcdefghijklmnopqrstuvwxyz'
CHARACTER(LEN=*),PARAMETER   :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
LOGICAL                      :: HasEq
!===================================================================================================================================
HasEq=.FALSE.
Str2=Str1
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  ! Transformation stops at "="
  IF(Str1(iLen:iLen).EQ.'=') HasEq=.TRUE.
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF ((Upper > 0).AND. .NOT. HasEq) Str2(iLen:iLen) = lc(Upper:Upper)
END DO
END SUBROUTINE LowCase

END MODULE MOD_ReadInTools
