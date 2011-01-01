MODULE randomga
    USE IFPORT
    USE parameters
    USE statistics
    IMPLICIT NONE

CONTAINS

    SUBROUTINE setSeed()
        IMPLICIT NONE
        INTEGER (KIND=integerSize)      :: allocationError
        INTEGER (KIND=4), ALLOCATABLE    :: values(:) !Contains values from the time/date call

        ALLOCATE(values(8), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("values")

        CALL DATE_AND_TIME(VALUES=values) !These values depend on the current time

        IF(values(8)==0) THEN                !We only want a non-zero value
             values(8)=values(5)+values(6)+values(7)+values(8) !Substitute in the case of zero (HH+MM+SS)
        END IF

        IF(randomSeed==-1) THEN
            CALL SRAND(values(8))
            randomSeed = values(8)
        ELSE
            CALL SRAND( INT(randomSeed) )
        END IF

        DEALLOCATE(values)
    END SUBROUTINE

    LOGICAL (KIND=logicalSize) FUNCTION getRandomLogical()
        IMPLICIT NONE
        getRandomLogical = MOD( IRAND(), 2 )                     !GENERATE RANDOM .TRUE. OR FALSE
    END FUNCTION

    INTEGER (KIND=integerSize) FUNCTION getIntRandomMax(maximumNumber)
        IMPLICIT NONE
        INTEGER (KIND=integerSize), INTENT(IN) :: maximumNumber
        getIntRandomMax = MOD( IRAND(), maximumNumber ) + 1                 !GENERATE RANDOM FROM 1 TO maximumNumber
    END FUNCTION

    INTEGER (KIND=integerSize) FUNCTION getIntRandomMinMax(minimumNumber, maximumNumber)
        IMPLICIT NONE
        INTEGER (KIND=integerSize), INTENT(IN) :: minimumNumber, maximumNumber
        getIntRandomMinMax = MOD( IRAND(), maximumNumber-minimumNumber+1 ) + minimumNumber     !GENERATE RANDOM FROM minimumNumber TO maximumNumber
    END FUNCTION

    REAL (KIND=realSize) FUNCTION getRealRandomMax(maximumNumber)
        IMPLICIT NONE
        REAL (KIND=realSize), INTENT(IN) :: maximumNumber
        getRealRandomMax = MOD( RAND(), maximumNumber )                     !GENERATE RANDOM FROM 0 TO maximumNumber
    END FUNCTION

    REAL (KIND=realSize) FUNCTION getRealRandom()
        IMPLICIT NONE
        getRealRandom = RAND()                                     !GENERATE RANDOM FROM 0 TO 1
    END FUNCTION

END MODULE randomga

