MODULE selectors
    USE parameters
    USE sorter
    USE randomga
    IMPLICIT NONE !Não permite variaves não declaradas


CONTAINS

!-------------------------------------------------------------------------------
    SUBROUTINE selectPopTTB
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i, counter
        counter = 1
        DO i=1, offspringSize
            IF(i>popKeep) counter = 1
            geneForMate(i) = counter
            counter = counter + 1
        END DO
    END SUBROUTINE


!-------------------------------------------------------------------------------


!RANKWEIGHTING------------------------------------------------------------------
    SUBROUTINE selectPopRankWeighting
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i
        DO i=1, offspringSize
            geneForMate(i) = selectGeneRankWeighting()
        END DO
    END SUBROUTINE

    FUNCTION selectGeneRankWeighting()
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)                 :: selectGeneRankWeighting
        INTEGER (KIND=integerSize), SAVE         :: lastMate = 0

        CALL selectMate(selectGeneRankWeighting)

        IF( selectGeneRankWeighting==lastMate ) selectGeneRankWeighting = getIntRandomMax(popSize)

        lastMate = selectGeneRankWeighting

    CONTAINS
        SUBROUTINE selectMate(mate)
            IMPLICIT NONE !Não permite variaves não declaradas
            INTEGER (KIND=integerSize), INTENT(OUT) :: mate
            INTEGER (KIND=integerSize)                 :: i
            REAL (KIND=realSize)                    :: tmp
            tmp = getRealRandom()

            IF( tmp>=0 .AND. tmp<rankweighProbability(1) ) THEN
                mate = 1
                RETURN
            END IF

            DO i=1, popKeep-1
                IF( tmp>=rankweighProbability(i) .AND. tmp<rankweighProbability(i+1) ) THEN
                    mate = i+1
                    RETURN
                END IF
            END DO
        END SUBROUTINE
    END FUNCTION

    SUBROUTINE calculateRankweigh
        INTEGER (KIND=integerSize) :: i
        REAL (KIND=realSize) :: rankweighProbabilityTmp( popKeep )
        REAL (KIND=realSize) :: soma = 0.0
        DO i=1, popKeep
            soma = soma + REAL(i)
        END DO
        DO i=1, popKeep
            rankweighProbabilityTmp(i) = REAL((popKeep - i + 1)/REAL(soma))
        END DO
        DO i=1, popKeep
            rankweighProbability(i) = SUM(rankweighProbabilityTmp(1:i))
        END DO
    END SUBROUTINE
!-------------------------------------------------------------------------------


!selectPopRandom----------------------------------------------------------------
    SUBROUTINE selectPopRandom
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i

        DO i=1, offspringSize
            geneForMate(i) = selectGeneRandom()
        END DO
    END SUBROUTINE

    FUNCTION selectGeneRandom()
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)                 :: selectGeneRandom
        INTEGER (KIND=integerSize), SAVE         :: lastMate

        selectGeneRandom = getIntRandomMax(popKeep)

        IF( selectGeneRandom==lastMate ) selectGeneRandom = getIntRandomMax(popSize)

        lastMate = selectGeneRandom

    END FUNCTION
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE selectPopTournament
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i
        DO i=1, offspringSize
            geneForMate(i) = selectGeneTournament()
        END DO
    END SUBROUTINE

    FUNCTION selectGeneTournament()
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)    :: selectGeneTournament
        INTEGER (KIND=integerSize)     :: subset(2)

        subset(1) = getIntRandomMax(popSize)
        subset(2) = getIntRandomMax(popSize)

        IF(getRealRandom() <= tournamentSelectionRate) THEN
            IF(geneScore(subset(1)) > geneScore(subset(2))) THEN
                selectGeneTournament = subset(1)
            ELSE
                selectGeneTournament = subset(2)
            END IF
        ELSE
            IF(geneScore(subset(1)) < geneScore(subset(2))) THEN
                selectGeneTournament = subset(1)
            ELSE
                selectGeneTournament = subset(2)
            END IF
        END IF
    END FUNCTION
!-------------------------------------------------------------------------------


!CostWeighting------------------------------------------------------------------
    SUBROUTINE selectPopCostWeighting
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i
        CALL calcCostWeighting
        DO i=1, offspringSize
            geneForMate(i) = selectGeneCostWeighting()
        END DO
    CONTAINS
        SUBROUTINE calcCostWeighting
            IMPLICIT NONE !Não permite variaves não declaradas
            INTEGER (KIND=integerSize)     :: i
            REAL (KIND=realSize)        :: costWeighProbabilityTmp( popKeep )
            REAL (KIND=realSize)         :: soma
            soma = SUM( geneScore(1:popKeep) )
            DO i=1, popKeep
                costWeighProbabilityTmp(i) = geneScore(i)/soma
                costWeighProbability(i) = SUM(costWeighProbabilityTmp(1:i))
            END DO
        END SUBROUTINE
    END SUBROUTINE

    FUNCTION selectGeneCostWeighting()
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)                 :: selectGeneCostWeighting
        INTEGER (KIND=integerSize), SAVE         :: lastMate = 0

        CALL selectMate(selectGeneCostWeighting)

        IF( selectGeneCostWeighting==lastMate ) selectGeneCostWeighting = getIntRandomMax(popSize)

        lastMate = selectGeneCostWeighting

    CONTAINS
        SUBROUTINE selectMate(mate)
            IMPLICIT NONE !Não permite variaves não declaradas
            INTEGER (KIND=integerSize), INTENT(OUT) :: mate
            INTEGER (KIND=integerSize)                 :: i
            REAL (KIND=realSize)                    :: tmp
            tmp = getRealRandom()

            IF( tmp>=0 .AND. tmp<costWeighProbability(1) ) THEN
                mate = 1
                RETURN
            END IF

            DO i=1, popKeep-1
                IF( tmp>=costWeighProbability(i) .AND. tmp<costWeighProbability(i+1) ) THEN
                    mate = i+1
                    RETURN
                END IF
            END DO
        END SUBROUTINE
    END FUNCTION
!-------------------------------------------------------------------------------




END MODULE
