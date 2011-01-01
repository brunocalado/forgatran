MODULE cross
    USE parameters
    USE randomga
    IMPLICIT NONE

CONTAINS

    !> Cross two genes in a single point\n\n
    !! Cross point: 3\n
    !! O O O | O O O <- parent 1\n
    !! X X X | X X X <- parent 2\n
    !! O O O | X X X <- offspring 1\n
    !! X X X | O O O <- offspring 2\n\n
    !! If niche is enable it will split the gene in n pieces
    !! and make single point cross piece by piece
    !! || -> Represents the niche size
    !! | -> Represents the cross point
    !! O | O O || O O | O || O | O O <- parent 1\n
    !! X | X X || X X | O || X | X X <- parent 2\n
    !! O | O O || O O | X || O | X X <- offspring 1\n
    !! O | O X || X X | O || O | O O <- offspring 2\n\n
    SUBROUTINE singleCross
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i
        INTEGER (KIND=integerSize)  :: crossPoint

        DO i=1, offspringSize/2, 2
            mask  = .FALSE.
            crossPoint = getIntRandomMax(totalGeneSize-1) !GENERATE RANDOM FROM 1 TO geneSize-1
            mask(1:crossPoint) = .TRUE.

            popAux(:, i)   = MERGE(pop(:, geneForMate(i) ), pop(:, geneForMate(i+1)), mask)
            popAux(:, i+1) = MERGE(pop(:, geneForMate(i+1) ), pop(:, geneForMate(i)), mask)
        END DO

    END SUBROUTINE

    !> Cross two genes in a single point\n\n
    !! Cross point: 3\n
    !! O O O | O O O <- parent 1\n
    !! X X X | X X X <- parent 2\n
    !! O O O | X X X <- offspring 1\n
    !! X X X | O O O <- offspring 2\n\n
    !! If niche is enable it will split the gene in n pieces
    !! and make single point cross piece by piece
    !! || -> Represents the niche size
    !! | -> Represents the cross point
    !! O | O O || O O | O || O | O O <- parent 1\n
    !! X | X X || X X | O || X | X X <- parent 2\n
    !! O | O O || O O | X || O | X X <- offspring 1\n
    !! O | O X || X X | O || O | O O <- offspring 2\n\n
    SUBROUTINE singleCrossNiche
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i, j
        INTEGER (KIND=integerSize)  :: nichePointA, nichePointB
        INTEGER (KIND=integerSize)  :: crossPoint

        DO i=1, offspringSize/2, 2
            DO j=1, yGeneSize
                maskNiche  = .FALSE.
                crossPoint = getIntRandomMax(geneSize-1) !GENERATE RANDOM FROM 1 TO bitsForEncodedReal-1
                maskNiche(1:crossPoint) = .TRUE.

                nichePointA = geneSize*(j-1)+1
                nichePointB = nichePointA+geneSize-1

                popAux(nichePointA:nichePointB, i) = MERGE(                         &
                pop(nichePointA:nichePointB, geneForMate(i) ),                      &
                pop(nichePointA:nichePointB, geneForMate(i+1)), maskNiche)

                popAux(nichePointA:nichePointB, i+1) = MERGE(                       &
                pop(nichePointA:nichePointB, geneForMate(i+1) ),                    &
                pop(nichePointA:nichePointB, geneForMate(i)), maskNiche)

            END DO
        END DO

    END SUBROUTINE

    !> Cross two genes in two points\n\n
    !! Sample
    !! Cross points: 4 and 6\n
    !! O O O O | O O | O O O <- parent 1\n
    !! X X X X | X X | X X X <- parent 2\n
    !! O O O O | X X | O O O <- offspring 1\n
    !! X X X X | O O | X X X <- offspring 2
    SUBROUTINE doubleCross
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i
        INTEGER (KIND=integerSize)  :: crossPointA, crossPointB

        DO i=1, offspringSize/2, 2
            mask  = .FALSE.

            crossPointA = getIntRandomMinMax(2, totalGeneSize-2)
            crossPointB = getIntRandomMinMax(crossPointA+1, totalGeneSize-1)

            mask(crossPointA:crossPointB) = .TRUE.

            popAux(:, i)    = MERGE(pop(:, geneForMate(i) ), pop(:, geneForMate(i+1)), mask)
            popAux(:, i+1)  = MERGE(pop(:, geneForMate(i+1) ), pop(:, geneForMate(i)), mask)
        END DO

    END SUBROUTINE

    !> Cross two genes in two points\n\n
    !! Sample
    !! Cross points: 4 and 6\n
    !! O O O O | O O | O O O <- parent 1\n
    !! X X X X | X X | X X X <- parent 2\n
    !! O O O O | X X | O O O <- offspring 1\n
    !! X X X X | O O | X X X <- offspring 2
    SUBROUTINE doubleCrossNiche
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i, j
        INTEGER (KIND=integerSize)  :: crossPointA, crossPointB
        INTEGER (KIND=integerSize)  :: nichePointA, nichePointB

        DO i=1, offspringSize/2, 2
            DO j=1, yGeneSize
                maskNiche  = .FALSE.
                crossPointA = getIntRandomMinMax(2, geneSize-2)
                crossPointB = getIntRandomMinMax(crossPointA+1, geneSize-1)

                maskNiche(crossPointA:crossPointB) = .TRUE.

                nichePointA = geneSize*(j-1)+1
                nichePointB = nichePointA+geneSize-1

                popAux(nichePointA:nichePointB, i) = MERGE(                         &
                pop(nichePointA:nichePointB, geneForMate(i) ),                      &
                pop(nichePointA:nichePointB, geneForMate(i+1)), maskNiche)

                popAux(nichePointA:nichePointB, i+1) = MERGE(                       &
                pop(nichePointA:nichePointB, geneForMate(i+1) ),                    &
                pop(nichePointA:nichePointB, geneForMate(i)), maskNiche)

            END DO
        END DO
    END SUBROUTINE


    !> Cross two genes using Clever TSP\n\n
    !! 1 2 3 | 4 5 6 | 7 8 9 <- parent 1\n
    !! X X X | X X X | X X X <- parent 2\n
    !! O O O | X X X | O O O <- offspring 1\n
    !! X X X | O O O | X X X <- offspring 2
    SUBROUTINE permutationCross
        IMPLICIT NONE !Não permite variaves não declaradas
        LOGICAL (KIND=logicalSize) :: flag
        INTEGER (KIND=integerSize) :: i, j, k, l 
        INTEGER (KIND=integerSize) :: counter
        INTEGER (KIND=integerSize) :: crossPointA, crossPointB

        DO i=1, offspringSize/2, 2
    
            DO l=0, 1                                  

!                do delta=1, 10
!                    pop(delta,geneForMate(i+l)) = delta
!                end do
!                do delta=10, 1, -1
!                    pop(delta,geneForMate(i+1-l)) =11- delta
!                end do
!
!                print *, pop(:,geneForMate(i+l))
!                print *, pop(:,geneForMate(i+1-l))


                crossPointA = getIntRandomMinMax(1, geneSize-2)
                crossPointB = getIntRandomMinMax(crossPointA+1, geneSize)

                counter = crossPointB - crossPointA + 1 

                popAux(1:counter, i+l) = pop(crossPointA:crossPointB, geneForMate(i+l))

!                print *
!                print *,  popAux(1:counter, i+l)
!                print *

                counter = counter + 1
                DO j=1, totalGeneSize
                    flag = .TRUE.
                    DO k=crossPointA, crossPointB
                        IF(pop(j, geneForMate(i+1-l)) == pop(k, geneForMate(i+l))) THEN
                            flag = .FALSE.
                            EXIT
                        END IF
                    END DO

                    IF(flag) THEN
                        popAux(counter, i+l) = pop(j, geneForMate(i+1-l))
                        counter = counter + 1
                    END IF
                END DO

!                print *, popAux(: , i+l)  
!
!                print *
!                print *
!                print *
!
            END DO

!
!                pause



        END DO


    END SUBROUTINE


    !> Cross two genes in two points in the block border\n\n
    !! This will choose just crosspoints blockSize multipliers\n\n
    !! Sample\n
    !! blockSize = 3\n
    !! Cross points: 3 and 6\n
    !! O O O | O O O | O O O <- parent 1\n
    !! X X X | X X X | X X X <- parent 2\n
    !! O O O | X X X | O O O <- offspring 1\n
    !! X X X | O O O | X X X <- offspring 2
    SUBROUTINE doubleCrossPermutation(blockSize)
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize), INTENT(IN)  :: blockSize
        INTEGER (KIND=integerSize)              :: crossPointA, crossPointB
        INTEGER (KIND=integerSize)              :: i, totalBlocks
        totalBlocks = geneSize/blockSize

        DO i=1, offspringSize/2, 2
            mask  = .FALSE.

            crossPointA = getIntRandomMinMax(2, totalBlocks-2)
            crossPointB = getIntRandomMinMax(crossPointA+1, totalBlocks-1)

            crossPointA = 1 + blockSize*(crossPointA-1)
            crossPointB = blockSize + blockSize*(crossPointB-1)

            mask(crossPointA:crossPointB) = .TRUE.

            popAux(:, i)    = MERGE(pop(:, geneForMate(i) ), pop(:, geneForMate(i+1)), mask)
            popAux(:, i+1)  = MERGE(pop(:, geneForMate(i+1) ), pop(:, geneForMate(i)), mask)
        END DO
    END SUBROUTINE


    !> Cross two genes in poliCrossSize points\n\n 
    !! Sample\n
    !! poliCrossSize = 3\n
    !! Cross points: 2, 4 and 6\n
    !! O O | O O | O O | O O O <- parent 1\n
    !! X X | X X | X X | X X X <- parent 2\n
    !! O O | X X | O O | X X X <- offspring 1\n
    !! X X | O O | X X | O O O <- offspring 2
    SUBROUTINE poliCross
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i, j
        INTEGER (KIND=integerSize)  :: crossPoint(poliCrossSize)

        DO i=1, offspringSize/2, 2
            mask = .FALSE.
            crossPoint(1) = getIntRandomMax(geneSize-poliCrossSize)
            mask(1:crossPoint(1)) = .TRUE.
            DO j=2, poliCrossSize
                crossPoint(j) = getIntRandomMinMax(crossPoint(j-1)+1, geneSize-(poliCrossSize+1-j))
            END DO

            DO j=2, poliCrossSize, 2
                IF(j==poliCrossSize) THEN
                    mask(crossPoint(j):) = .TRUE.
                    EXIT
                END IF
                mask(crossPoint(j):crossPoint(j+1)) = .TRUE.
            END DO

            popAux(:, i)    = MERGE(pop(:, geneForMate(i) ), pop(:, geneForMate(i+1)), mask)
            popAux(:, i+1)  = MERGE(pop(:, geneForMate(i+1) ), pop(:, geneForMate(i)), mask)
        END DO
    END SUBROUTINE

    !> Cross two genes using uniformCross
    SUBROUTINE uniformCross
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i, j
        DO i=1, offspringSize/2, 2
            mask = .FALSE.
            mask = (/(getRandomLogical() ,j=1, totalGeneSize)/)
            popAux(:, i) = MERGE(pop(:, geneForMate(i) ), pop(:, geneForMate(i+1)), mask)
            popAux(:, i+1) = MERGE(pop(:, geneForMate(i+1) ), pop(:, geneForMate(i)), mask)
        END DO

    END SUBROUTINE

    !test me
    SUBROUTINE uniformCrossNiche
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i, j, k
        INTEGER (KIND=integerSize)  :: nichePointA, nichePointB
        INTEGER (KIND=integerSize)  :: crossPoint

        DO i=1, offspringSize/2, 2
            DO j=1, yGeneSize
                maskNiche  = .FALSE.
                maskNiche = (/(getRandomLogical() ,k=1, yGeneSize)/)
                
                nichePointA = geneSize*(j-1)+1
                nichePointB = nichePointA+geneSize-1

                popAux(nichePointA:nichePointB, i) = MERGE(                         &
                pop(nichePointA:nichePointB, geneForMate(i) ),                      &
                pop(nichePointA:nichePointB, geneForMate(i+1)), maskNiche)

                popAux(nichePointA:nichePointB, i+1) = MERGE(                       &
                pop(nichePointA:nichePointB, geneForMate(i+1) ),                    &
                pop(nichePointA:nichePointB, geneForMate(i)), maskNiche)
    
            END DO
        END DO

    END SUBROUTINE

    !> Cross two genes using uniformCross
    SUBROUTINE uniformCrossPlus
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize) :: i, j
        INTEGER (KIND=integerSize) :: trueBitsPerGene
        trueBitsPerGene = CEILING(0.6*totalGeneSize)

        DO i=1, offspringSize/2, 2
            mask = .FALSE.
            DO j=1, trueBitsPerGene
                mask(getIntRandomMax(totalGeneSize)) = .TRUE.
            END DO
            popAux(:, i) = MERGE(pop(:, geneForMate(i) ), pop(:, geneForMate(i+1)), mask)
            popAux(:, i+1) = MERGE(pop(:, geneForMate(i+1) ), pop(:, geneForMate(i)), mask)

        END DO
    END SUBROUTINE

    !test me
    !> Cross two genes using uniformCross
    SUBROUTINE uniformCrossPlusNiche
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER (KIND=integerSize)  :: i, j, k
        INTEGER (KIND=integerSize)  :: nichePointA, nichePointB
        INTEGER (KIND=integerSize) :: trueBitsPerGene
        trueBitsPerGene = CEILING(0.6*totalGeneSize)

        DO i=1, offspringSize/2, 2
            DO j=1, yGeneSize
                maskNiche  = .FALSE.
                DO k=1, trueBitsPerGene
                    mask(getIntRandomMax(yGeneSize)) = .TRUE.
                END DO
                
                nichePointA = geneSize*(j-1)+1
                nichePointB = nichePointA+geneSize-1

                popAux(nichePointA:nichePointB, i) = MERGE(                         &
                pop(nichePointA:nichePointB, geneForMate(i) ),                      &
                pop(nichePointA:nichePointB, geneForMate(i+1)), maskNiche)

                popAux(nichePointA:nichePointB, i+1) = MERGE(                       &
                pop(nichePointA:nichePointB, geneForMate(i+1) ),                    &
                pop(nichePointA:nichePointB, geneForMate(i)), maskNiche)

            END DO
        END DO

    END SUBROUTINE

END MODULE
