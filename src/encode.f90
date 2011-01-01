MODULE encode
    USE parameters
    IMPLICIT NONE !Não permite variaves não declaradas


CONTAINS
    INTEGER (KIND=integerSize) FUNCTION binaryToInteger(gene, parameterPosition)
        LOGICAL (KIND=logicalSize), INTENT(IN)     :: gene(:)
        INTEGER (KIND=integerSize), INTENT(IN)     :: parameterPosition
        INTEGER (KIND=integerSize)                 :: i
        INTEGER (KIND=integerSize)                 :: nichePointA
        nichePointA = geneSize*(parameterPosition-1)+1

        binaryToInteger = 0
        DO i=0, geneSize-1
            IF(gene(nichePointA+i))    binaryToInteger = binaryToInteger + 2**i
        END DO

    END FUNCTION

    
    !FIX ME
    FUNCTION integerToBinary(value, length)
        INTEGER (KIND=integerSize), INTENT(IN)     :: value
        INTEGER (KIND=integerSize), INTENT(IN)     :: length
        LOGICAL (KIND=logicalSize)                 :: integerToBinary(length)
        INTEGER (KIND=integerSize)                 :: i, tmp
        integerToBinary = .FALSE.

        IF(value==0) RETURN
        IF(value==1) THEN
            integerToBinary(1) = .TRUE.
            RETURN
        END IF

        tmp=value
        DO i=1, length
            IF(MOD(tmp, 2)==1) integerToBinary(i) = .TRUE.

            IF(tmp/2==1) THEN
                integerToBinary(i+1) = .TRUE.
                RETURN
            END IF
            tmp = tmp / 2
        END DO
    END FUNCTION

    REAL (KIND=realSize) FUNCTION binaryToReal(gene, parameterPosition)
        IMPLICIT NONE !Não permite variaves não declaradas
        LOGICAL (KIND=logicalSize), INTENT(IN)     :: gene(:)
        INTEGER (KIND=integerSize), INTENT(IN)     :: parameterPosition

        binaryToReal = lowerBoundForEncodedReal(parameterPosition) + &
        incrementForEncodedReal(parameterPosition)*REAL(binaryToInteger(gene, parameterPosition))
    END FUNCTION

    FUNCTION realToBinary(value, parameterPosition)
        REAL (KIND=realSize), INTENT(IN)         :: value
        INTEGER (KIND=integerSize), INTENT(IN)   :: parameterPosition
        LOGICAL (KIND=logicalSize)               :: realToBinary(geneSize)
        realToBinary=.TRUE.
    END FUNCTION
END MODULE
