MODULE converge
    USE parameters
    USE selectors
    IMPLICIT NONE
    INTEGER (KIND=integerSize), PARAMETER     :: setSize = 100
    INTEGER (KIND=integerSize), PARAMETER     :: deviationSize = 5
    REAL (KIND=realSize)                    :: scorePerBestGene(setSize)
    REAL (KIND=realSize)                    :: deviationPerSet(deviationSize)
CONTAINS

    !>RETURNS TRUE IF IT CONVERGED
    LOGICAL (KIND=logicalSize) FUNCTION convergeCheck(bestGeneScore)
        REAL (KIND=realSize), INTENT(IN)      :: bestGeneScore
        INTEGER (KIND=integerSize)            :: i
        INTEGER (KIND=integerSize), SAVE    :: counter = 1
        INTEGER (KIND=integerSize), SAVE    :: counterDeviation = 1
        INTEGER (KIND=integerSize)             :: counterTmp = 0
        REAL (KIND=realSize)                :: mean=0
        REAL (KIND=realSize)                :: deviation=0

        scorePerBestGene(counter) = bestGeneScore
        counter = counter + 1

        IF(counter > setSize) THEN
            counter = 1
            DO i=1, setSize
                mean = mean + scorePerBestGene(i)
            END DO
            mean = mean/setSize
            DO i=1, setSize
                deviation =  deviation + (scorePerBestGene(i) - mean) * (scorePerBestGene(i) - mean)
            END DO
            deviation = SQRT(deviation/setSize)

            deviationPerSet(counterDeviation) = deviation
            counterDeviation = counterDeviation + 1
            IF(counterDeviation >= deviationSize ) THEN
                counterDeviation = 1
                DO i=1, deviationSize
                    IF(deviationPerSet(i) == 0) THEN
                        counterTmp = counterTmp + 1
                    END IF
                END DO
                IF(counterTmp == deviationSize) THEN
                    convergeCheck = .TRUE.
                    RETURN
                END IF
            END IF
        END IF

        convergeCheck = .FALSE.
    END FUNCTION

    !>RETURNS TRUE IF IT CONVERGED
    LOGICAL (KIND=logicalSize) FUNCTION convergeCheckByThreshold( threshold )
        REAL (KIND=realSize), INTENT(IN) :: threshold

        IF( MAXVAL(geneScore) >= threshold) THEN
            convergeCheckByThreshold = .TRUE.
        ELSE
            convergeCheckByThreshold = .FALSE.
        END IF
    END FUNCTION
END MODULE


