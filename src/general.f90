MODULE general
    USE IFPORT
    USE parameters
    USE selectors
    USE score
    USE sorter
    USE randomga

!OPTIONAL----------------------------------------------------------------------
    USE statistics
    IMPLICIT NONE

CONTAINS
    !> Score the population from n to popSize
    SUBROUTINE scorePop(startPoint, endPoint)
        INTEGER (KIND=integerSize), INTENT(IN)  :: startPoint, endPoint
        INTEGER (KIND=integerSize)              :: i
        DO i=startPoint, endPoint
            geneScore(i) = scoreGene(i)
        END DO
    END SUBROUTINE

    SUBROUTINE allocateVars
        IMPLICIT NONE
        INTEGER    (KIND=integerSize)            :: allocationError

        popKeep         = NINT(popSize*selectionRate)
        offspringSize   = popSize-popKeep
        totalGeneSize   = geneSize*yGeneSize

        ALLOCATE(geneForMate(popSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("geneForMate")

        ALLOCATE(rankweighProbability(popKeep), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("rankweighProbability")

        ALLOCATE(costWeighProbability(popKeep), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("costWeighProbability")

        ALLOCATE(geneScore(popSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("geneScore")

        ALLOCATE(pop(totalGeneSize, popSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("pop")

        ALLOCATE(popAux(totalGeneSize, popSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("popAux")

        ALLOCATE(mask(totalGeneSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("mask")

        ALLOCATE(maskNiche(geneSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("maskNiche")

    END SUBROUTINE

    SUBROUTINE deallocateVars
        DEALLOCATE(geneForMate)
        DEALLOCATE(rankweighProbability)
        DEALLOCATE(costWeighProbability)
        DEALLOCATE(geneScore)
        DEALLOCATE(pop)
        DEALLOCATE(popAux)
        DEALLOCATE(mask)
        DEALLOCATE(maskNiche)
    END SUBROUTINE

    !>
    SUBROUTINE calculateRealParameters(lowerBounds, upperBounds)
        IMPLICIT NONE
        REAL (KIND=realSize), INTENT(IN) :: lowerBounds(:)
        REAL (KIND=realSize), INTENT(IN) :: upperBounds(:)
        INTEGER (KIND=integerSize)       :: i
        INTEGER (KIND=integerSize)       :: allocationError

        ALLOCATE(upperBoundForEncodedReal(yGeneSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("upperBoundForEncodedReal")

        ALLOCATE(lowerBoundForEncodedReal(yGeneSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("lowerBoundForEncodedReal")

        ALLOCATE(incrementForEncodedReal(yGeneSize), STAT=allocationError)
        IF(allocationError /= 0) CALL errorMessage("incrementForEncodedReal")

        lowerBoundForEncodedReal = lowerBounds
        upperBoundForEncodedReal = upperBounds

        !Calculate bits
        DO i=0, integerSize*8-1
            IF( 2**i > resolutionForEncodedReal ) THEN
                resolutionForEncodedReal = 2**i-1
                geneSize = i
                EXIT
            END IF
        END DO

        DO i=1, yGeneSize
            incrementForEncodedReal(i) = (                                        &
            ABS(upperBoundForEncodedReal(i) - lowerBoundForEncodedReal(i))        &
            ) / REAL(resolutionForEncodedReal)

        END DO


    END SUBROUTINE

!    SUBROUTINE calculateIntegerParameters(lowerBounds, upperBounds)
!        IMPLICIT NONE
!        INTEGER (KIND=integerSize), INTENT(IN)    :: lowerBounds(:)
!        INTEGER (KIND=integerSize), INTENT(IN)    :: upperBounds(:)
!        INTEGER (KIND=integerSize)                :: i, j
!        INTEGER    (KIND=integerSize)                :: allocationError
!
!        ALLOCATE(upperBoundForEncodedInteger(yGeneSize), STAT=allocationError)
!        IF(allocationError /= 0) CALL errorMessage("upperBoundForEncodedInteger")
!
!        ALLOCATE(lowerBoundForEncodedInteger(yGeneSize), STAT=allocationError)
!        IF(allocationError /= 0) CALL errorMessage("lowerBoundForEncodedInteger")
!
!        lowerBoundForEncodedInteger = lowerBounds
!        upperBoundForEncodedInteger = upperBounds
!
!        !Calculate bits
!        DO i=1, yGeneSize
!            DO j=0, integerSize*8-1
!                IF( 2**j > (upperBoundForEncodedInteger(i) - lowerBoundForEncodedInteger(i)) ) THEN
!                    geneSize = j
!                    EXIT
!                END IF
!            END DO
!        END DO
!
!        ALLOCATE(upperBoundForEncodedReal(yGeneSize), STAT=allocationError)
!        IF(allocationError /= 0) CALL errorMessage("upperBoundForEncodedReal")
!
!        ALLOCATE(lowerBoundForEncodedReal(yGeneSize), STAT=allocationError)
!        IF(allocationError /= 0) CALL errorMessage("lowerBoundForEncodedReal")
!
!        ALLOCATE(incrementForEncodedReal(yGeneSize), STAT=allocationError)
!        IF(allocationError /= 0) CALL errorMessage("incrementForEncodedReal")
!
!        lowerBoundForEncodedReal = lowerBounds
!        upperBoundForEncodedReal = upperBounds
!        incrementForEncodedReal = 1.0D0
!    END SUBROUTINE

    SUBROUTINE saveParamaters(file)
        IMPLICIT NONE
        CHARACTER (LEN=*, KIND=characterSize) , INTENT(IN)    :: file
        INTEGER (KIND=integerSize)                             :: iostatus

        NAMELIST /gaparameters/ geneSize, yGeneSize, popSize, maxGenerations, selectionRate,     &
        mutateRate, tournamentSelectionRate, randomSeed, poliCrossSize, resolutionForEncodedReal


        OPEN( UNIT=11, FILE=file, ACTION="WRITE", IOSTAT=iostatus )

        WRITE (11, NML = gaparameters)

        CLOSE (11, STATUS = 'KEEP')
    END SUBROUTINE

    SUBROUTINE loadParamaters(file)
        IMPLICIT NONE
        CHARACTER (LEN=*, KIND=characterSize) , INTENT(IN)    :: file
        INTEGER (KIND=integerSize)                             :: iostatus

        NAMELIST /gaparameters/ geneSize, yGeneSize, popSize, maxGenerations, selectionRate,     &
        mutateRate, tournamentSelectionRate, randomSeed, poliCrossSize, resolutionForEncodedReal


        OPEN( UNIT=11, FILE=file, ACTION="READ", IOSTAT=iostatus )

        READ (11, NML = gaparameters)

        CLOSE (11, STATUS = 'KEEP')

        WRITE (*, gaparameters)
    END SUBROUTINE

    SUBROUTINE variableMutationByFitness
        IMPLICIT NONE
        REAL (KIND=realSize) :: maximumScore
        REAL (KIND=realSize) :: mean
        REAL (KIND=realSize) :: delta
        maximumScore = MAXVAL(geneScore(1:popSize))
        mean = SUM(geneScore(1:popSize))/popSize
        IF( ISNAN(mean) ) RETURN
        delta = ((maximumScore - mean)/(maximumScore + mean))
        IF( ISNAN(delta) ) RETURN

        !- - - DEC$ IF debug .EQ. 1
        !    PRINT "(4A10)", "delta - ", "mutateRate - ", "maximumScore - ", "mean - "
        !    PRINT "(4F10.4)", delta, mutateRate, maximumScore, mean
        ! -  - DEC$ ENDIF


        IF( ABS(delta) <= lowestDelta .AND. mutateRate < highestMutateRate ) THEN
            mutateRate = mutateRate*mutateRateMultiplier
        ELSE IF( delta >= highestDelta .AND. mutateRate > lowestMutateRate ) THEN
            mutateRate = mutateRate/mutateRateMultiplier
        END IF

    END SUBROUTINE

END MODULE
