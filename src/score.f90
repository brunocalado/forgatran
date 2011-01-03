MODULE score
!NEEDED------------------------------------------------------------------------
    USE parameters
!OPTIONAL----------------------------------------------------------------------
    IMPLICIT NONE !Não permite variaves não declaradas

    !Your variables and constants
    INTEGER, PRIVATE :: i

CONTAINS
	!Binary Sample
    REAL (KIND=realSize) FUNCTION scoreGene(gene)
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize), INTENT(IN) :: gene
        scoreGene = 0.0D0

        !Put your fitness code here

        DO i=1, totalGeneSize, 2
            IF( pop(i, gene) == .TRUE. .AND. (pop(i+1, gene) == .FALSE.) ) scoreGene = scoreGene + 1.0D0
        END DO

    END FUNCTION

    !Integer Sample
!    REAL (KIND=realSize) FUNCTION scoreGene(gene)
!        IMPLICIT NONE !Não permite variaves não declaradas
!        INTEGER (KIND=integerSize), INTENT(IN) :: gene
!        scoreGene = 0.0D0
!
!        !Put your fitness code here
!        DO i=1, totalGeneSize, 2
!            IF( pop(i, gene) == 1 .AND. (pop(i+1, gene) == 2) ) scoreGene = scoreGene + 1.0D0
!        END DO
!
!    END FUNCTION
 
END MODULE



