PROGRAM main
!SELECT THE GENE---------------------------------------------------------------
     USE binaryGene
!    USE integerGene
!NEEDED------------------------------------------------------------------------
    USE selectors
    USE cross
    USE replacement
    USE general
!OPTIONAL-----------------------------------------------------------------------
    USE converge
    USE score
    USE statistics
    USE encode
!
 !------------------------------------------------------------------------------
    IMPLICIT NONE !Não permite variaves não declaradas
    INTEGER (KIND=integerSize), PARAMETER    :: boundsSize = 2
    REAL (KIND=realSize), PARAMETER            :: upperBounds(boundsSize) = 5.0d0
    REAL (KIND=realSize), PARAMETER            :: lowerBounds(boundsSize) = -5.0d0

!    CALL loadParamaters("forgatran.config")

!    CALL calculateRealParameters(lowerBounds, upperBounds)           !Call this if you will use binary real encode


    CALL startup

    CALL createPop
!    CALL loadPop("pop.data")

!    CALL printHeader


    CALL startCounter
    CALL run
    CALL finishCounter


    CALL savePop("pop.data")

    CALL viewPop(04)
!    CALL show



    CALL shutdown
CONTAINS
    SUBROUTINE run
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i

        CALL scorePop(1, popSize)

        DO i=1, maxGenerations

!SORT---------------------------------------------------------------------------
                CALL sort(geneScore)


!MUTATION VARIATION-------------------------------------------------------------
                CALL variableMutationByFitness

!SELECT-------------------------------------------------------------------------
!                CALL selectPopTTB !TOP TO BOTTOM
!                CALL selectPopRandom
                CALL selectPopRankWeighting !Roulette Wheel by Rank Weighting
!                CALL selectPopCostWeighting
!                CALL selectPopTournament

!CROSS--------------------------------------------------------------------------
!                CALL singleCross
!                CALL singleCrossNiche
                CALL doubleCross
!                CALL doubleCrossNiche
!                CALL permutationCross
!                CALL doubleCrossPermutation(2)
!                CALL poliCross
!                CALL uniformCross
!                CALL uniformCrossNiche
!                CALL uniformCrossPlus
!                CALL uniformCrossPlusNiche


!MUTATE-------------------------------------------------------------------------
!                CALL mutate
                CALL mutatePlus

!REPLACEMENT-OR INSERTION-------------------------------------------------------
                CALL replaceWorst
!                CALL insert

!SCORE POP----------------------------------------------------------------------
                CALL scorePop(popKeep+1, popSize) !FROM popKeep+1 to popSize



!CONVERGENCE CHECK--------------------------------------------------------------
                IF( convergeCheckByThreshold( REAL(totalGeneSize)/2.0D0) ) THEN
                    PRINT "(A, I10)", "Convergence achieved at: ", i
                    RETURN
                END IF

!                IF(convergeCheck(geneScore(1))) THEN
!                    PRINT "(A, I10)", "Convergence achieved at: ", (i-1)*extendsGenerations+j
!                    RETURN
!                END IF


!OPTIONAL-----------------------------------------------------------------------
                CALL progress(i)
!                stat_scorePerGeneration(i) = averageScore(geneScore)
!                stat_scorePerBestGene(i) = geneScore(1)



!                CALL stat_mutate(i, geneScore)
        END DO

    END SUBROUTINE

    SUBROUTINE startup
        CALL allocateVars
        CALL setSeed
        CALL calculateRankweigh !RANKWEIGHTING - P(n) = (Nkeep - n + 1)/SUM(popKeep)
!        CALL allocateMean
    END SUBROUTINE

    SUBROUTINE shutdown
!        CALL deallocateMean
        CALL deallocateVars
    END SUBROUTINE
END PROGRAM

