MODULE replacement
    USE parameters
    USE general
    USE sorter
    IMPLICIT NONE !Não permite variaves não declaradas

CONTAINS
    SUBROUTINE replaceWorst
        IMPLICIT NONE
        pop(:, (popKeep+1):popSize) = popAux(:,1:offspringSize)
    END SUBROUTINE

    SUBROUTINE insert
    !TODO usar move_alloc para ficar mais rapido
        IMPLICIT NONE
        INTEGER (KIND=integerSize)              :: i
        LOGICAL (KIND=logicalSize)                :: popTmp(totalGeneSize, popSize)
        REAL (KIND=realSize)                    :: geneScoreTmp(popSize)

        popTmp(:,1:popSize) = pop(:,1:popSize)
        geneScoreTmp(1:popSize) = geneScore(1:popSize)

        popSize = popSize*2

        DEALLOCATE(geneScore)
        DEALLOCATE(pop)

        ALLOCATE(geneScore(popSize))
        ALLOCATE(pop(totalGeneSize, popSize))

        pop(:, 1:(popSize/2)) = popTmp(:,1:popSize)
        pop(:, (popSize/2+1):popSize) = popAux
        geneScore(1:(popSize/2)) = geneScoreTmp(1:popSize)

        CALL scorePop(popSize/2+1, popSize)

        CALL sort(geneScore)

        popSize = popSize/2

        popTmp = pop(:,1:popSize)
        geneScoreTmp = geneScore(1:popSize)

        DEALLOCATE(geneScore)
        DEALLOCATE(pop)

        ALLOCATE(geneScore(popSize))
        ALLOCATE(pop(totalGeneSize, popSize))

        pop(:, 1:popSize) = popTmp
        geneScore = geneScoreTmp

    END SUBROUTINE


END MODULE
