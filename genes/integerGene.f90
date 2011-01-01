!> Integer gene with one dimension. \n\n
!! Enable it with: USE integer in the main file.
MODULE integerGene
    USE parameters
    USE sorter
    USE randomga
    USE general

!OPTIONAL----------------------------------------------------------------------
    USE statistics
    IMPLICIT NONE

    !Auxiliary Variables


CONTAINS

    !> Create the population with random values.
    SUBROUTINE createPop
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i
        DO i=1, popSize
            pop(:, i) = createGene()
        END DO
    END SUBROUTINE

    FUNCTION createGene()
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)    :: createGene(totalGeneSize)
        INTEGER (KIND=integerSize)    :: i
        DO i=1, totalGeneSize
            createGene(i) = getIntRandomMinMax(minimumInteger, maximumInteger)
        END DO
    END FUNCTION

    !> Test every bit from the new offspring
    SUBROUTINE mutate
        IMPLICIT NONE
        INTEGER (KIND=integerSize)    :: i, j
        DO i=1, offspringSize
            DO j=1, totalGeneSize
                IF( getRealRandom() <= mutateRate ) THEN
                    popAux(j, i) = getIntRandomMinMax(minimumInteger, maximumInteger)
                END IF
            END DO
        END DO
    END SUBROUTINE

    SUBROUTINE mutatePlus
        IMPLICIT NONE
        INTEGER (KIND=integerSize)    :: i, j
        INTEGER (KIND=integerSize)    :: mutationsPerGene, mutationPos
        mutationsPerGene = CEILING(mutateRate*totalGeneSize)
        DO i=1, offspringSize
            DO j=1, mutationsPerGene
                mutationPos = getIntRandomMax(totalGeneSize)
                popAux(mutationPos, i) = getIntRandomMinMax(minimumInteger, maximumInteger)
            END DO
        END DO
    END SUBROUTINE



!------------------------------------------------------------------------------
    SUBROUTINE viewPop(n)
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize), INTENT(IN)    :: n
        INTEGER (KIND=integerSize)                :: i

        CALL scorePop(1, popSize)
        CALL sort(geneScore)

        PRINT *
        PRINT "(A)","Rank - Score"
        DO i=1, n
            PRINT "(I4,A5,F20.4)", i, "   - ", geneScore(i)
            PRINT 20, pop(:, i)
        END DO

        20 FORMAT( <geneSize*yGeneSize> I )
    END SUBROUTINE

    SUBROUTINE savePop(file)
        IMPLICIT NONE
        CHARACTER (LEN=*, KIND=characterSize), INTENT(IN)    :: file
        INTEGER (KIND=integerSize)                             :: i, iostatus

        CALL scorePop(1, popSize)
        CALL sort(geneScore)

        OPEN( UNIT=11, FILE=file, ACTION='WRITE', IOSTAT=iostatus )
        DO i=1, popSize
            IF(iostatus>0) THEN
                PRINT "(3A)", "Write: ", file, " file fail."
                RETURN
            END IF
            WRITE(11, 20) pop(:, i)
        END DO

        20 FORMAT( <geneSize*yGeneSize> I )

        CLOSE (11, STATUS = 'KEEP')
    END SUBROUTINE

    SUBROUTINE loadPop(file)
        IMPLICIT NONE
        CHARACTER (LEN=*, KIND=characterSize), INTENT(IN)    :: file
        INTEGER (KIND=integerGeneSize)                            :: gene(geneSize*yGeneSize)
        INTEGER (KIND=integerSize)                             :: i, iostatus

        OPEN( UNIT=11, FILE=file, ACTION='READ', IOSTAT=iostatus )
        DO i=1, popSize
            IF(iostatus>0) THEN
                PRINT "(3A)", "Read: ", file, " file fail. Random new population created."
                CALL createPop
                RETURN
            END IF

            READ(11, 20) gene
            pop(:, i) = gene
        END DO

        20 FORMAT( <geneSize*yGeneSize> I )

        CLOSE (11, STATUS = 'KEEP')
        CALL scorePop(1, popSize) ! FROM 1 to popSize
    END SUBROUTINE
!------------------------------------------------------------------------------




END MODULE

