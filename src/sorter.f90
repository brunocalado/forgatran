MODULE sorter
    USE parameters
    IMPLICIT NONE !Não permite variaves não declaradas

CONTAINS
    SUBROUTINE sort(geneScore)
        IMPLICIT NONE
        REAL (KIND=realSize), INTENT(INOUT)       :: geneScore(popSize)
!        LOGICAL (KIND=logicalSize), INTENT(INOUT) :: pop(totalGeneSize, popSize)
        INTEGER (KIND=integerSize)                :: i
        INTEGER (KIND=integerSize)                :: indexScore(popSize)
        indexScore = (/ (i, i=1, popSize) /)
        CALL quickSort(geneScore, indexScore, popSize)
        pop(:,:) = pop(:,indexScore)
    END SUBROUTINE

    RECURSIVE SUBROUTINE quickSort(A, indexScore, Asize)
        IMPLICIT NONE
        INTEGER (KIND=integerSize), INTENT(IN)         :: Asize
        INTEGER (KIND=integerSize), INTENT(INOUT)     :: indexScore(Asize)
        REAL (KIND=realSize), INTENT(INOUT)         :: A(Asize)

        INTEGER (KIND=integerSize)                    :: index_tmp !FOR swap
        INTEGER (KIND=integerSize)                     :: iq
        REAL (KIND=realSize)                        :: tmp !FOR swap

        IF( Asize>1 ) THEN
            CALL partition( iq )
            CALL quickSort( A(:iq-1), indexScore(:iq-1), iq-1 )
            CALL quickSort( A(iq:),    indexScore(iq:), Asize-iq+1 )
        END IF

    CONTAINS
        SUBROUTINE partition(marker)
            IMPLICIT NONE
            INTEGER (KIND=integerSize), INTENT(OUT)        :: marker
            INTEGER (KIND=integerSize)                    :: i, j
            REAL (KIND=realSize)                         :: x
            x = A(1)
            i = 0
            j = Asize + 1

            DO
                j = j - 1
                DO
                    IF( A(j)>=x ) EXIT !CHANGE TO <
                    j = j - 1
                END DO
                i = i + 1
                DO
                    IF( A(i)<=x ) EXIT !CHANGE TO >
                    i = i + 1
                END DO
                IF( i<j ) THEN
                    CALL swap(i,j)
                ELSEIF( i==j ) THEN
                    marker = i + 1
                    RETURN
                ELSE
                    marker = i
                    RETURN
                END IF
            END DO

        END SUBROUTINE

        SUBROUTINE swap(posA, posB)
            IMPLICIT NONE
            INTEGER (KIND=integerSize), INTENT(IN)    :: posA, posB
            tmp = A(posA) ; A(posA) = A(posB) ; A(posB) = tmp
            index_tmp = indexScore(posA) ; indexScore(posA) = indexScore(posB) ; indexScore(posB) = index_tmp
        END SUBROUTINE
    END SUBROUTINE

END MODULE


