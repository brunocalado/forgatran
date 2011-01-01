MODULE statistics
    USE IFPORT
    USE parameters

    IMPLICIT NONE !Não permite variaves não declaradas
    INTEGER (KIND=8),PARAMETER, PRIVATE    :: stat_integerSize    = 8
    INTEGER (KIND=8),PARAMETER, PRIVATE    :: stat_realSize    = 8

    !CROSSING------------------------------------------------------------------
    INTEGER (KIND=stat_integerSize)     :: stat_CrossPerGeneration = 0

    !GENERATIONS---------------------------------------------------------------
    INTEGER (KIND=stat_integerSize)     :: stat_generation

    !MUTATION -----------------------------------------------------------------
    INTEGER (KIND=stat_integerSize)     :: stat_maxMutationPerGeneration = 0
    INTEGER (KIND=stat_integerSize)     :: stat_minMutationPerGeneration = 1000000
    INTEGER (KIND=stat_integerSize)     :: stat_mutationsPerGeneration = 0
    REAL (KIND=stat_realSize)             :: stat_averageMutations = 0

    !SCORE---------------------------------------------------------------------
    REAL (KIND=stat_realSize), ALLOCATABLE             :: stat_scorePerGeneration(:)
    REAL (KIND=stat_realSize), ALLOCATABLE             :: stat_scorePerBestGene(:)

CONTAINS
    SUBROUTINE logStatus(message)
        IMPLICIT NONE
        CHARACTER (LEN=*), INTENT(IN)       :: message
        INTEGER (KIND=integerSize)          :: iostatus
        INTEGER (KIND=integerSize), SAVE    :: counter = 0
        LOGICAL (KIND=logicalSize), SAVE    :: first = .TRUE.


        IF(first) THEN
            iostatus = DELFILESQQ("log")
            first = .FALSE.
        END IF

        OPEN( UNIT=20, FILE="log", ACTION='WRITE', IOSTAT=iostatus )
        IF(iostatus>0) THEN
            PRINT "(A)", "Write file fail."
            RETURN
        END IF

        WRITE(20, "(I3, A3, A)") counter, " - ", message

        counter = counter + 1

        CLOSE (20, STATUS = 'KEEP')
    END SUBROUTINE

    SUBROUTINE errorMessage(name)
        IMPLICIT NONE
        CHARACTER (LEN=*, KIND=characterSize) , INTENT(IN)  :: name

        PRINT *, name, " CANNOT be allocated. Too big?"
        PRINT *, "Terminating program execution"
        CALL logStatus(name // " CANNOT be allocated")
        CALL EXIT(-1)
    END SUBROUTINE
    REAL (KIND=realSize) FUNCTION taskTime(flag)
        IMPLICIT NONE
        LOGICAL (KIND=logicalSize), INTENT(IN)    :: flag
        REAL (KIND=realSize), SAVE                :: time_begin, time_end

        IF(flag) THEN
            time_begin = 0
            time_end = 0
            CALL CPU_TIME ( time_begin )
        ELSE
            CALL CPU_TIME ( time_end )
            taskTime = time_end - time_begin !seconds
            time_begin = 0
            time_end = 0
        ENDIF
    END FUNCTION

    SUBROUTINE startCounter
        REAL (KIND=realSize) :: dummy
        dummy = taskTime(.TRUE.)
    END SUBROUTINE

    SUBROUTINE finishCounter
        PRINT "(A,F15.3,A10)", 'Total time: ', taskTime(.FALSE.), " [seconds]"
    END SUBROUTINE



    SUBROUTINE progress(i)
        INTEGER (KIND=integerSize), INTENT(IN)    :: i

        IF( MOD(i, CEILING(maxGenerations*0.01) ) == 0 ) THEN
            OPEN(unit=6,carriagecontrol='fortran')
            WRITE(6, '($,A,A,F5.1,A6)') '+', CHAR(13), (REAL(i)/REAL(maxGenerations)) * 100, "%/100%"

!            PRINT "(F5.1,A6)", (REAL(i)/REAL(maxGenerations)) * 100,"%/100%"
        END IF


    END SUBROUTINE

    REAL (KIND=realSize) FUNCTION averageScore(score)
        IMPLICIT NONE
        REAL (KIND=realSize), INTENT(IN)    :: score(popSize)
        INTEGER (KIND=stat_integerSize)     :: i
        REAL (KIND=stat_realSize)             :: average
        average = 0
        DO i=1, popSize
            average = average + score(i)
        END DO
        averageScore = average/popSize
    END FUNCTION

    SUBROUTINE allocateMean
        ALLOCATE(stat_scorePerGeneration(maxGenerations), stat_scorePerBestGene(maxGenerations))
    END SUBROUTINE

    SUBROUTINE deallocateMean
!        INTEGER (KIND=integerSize) :: i, iostatus

!        OPEN( UNIT=11, FILE='stat_scorePerGeneration.data', ACTION='WRITE', IOSTAT=iostatus )
!        DO i=1, maxGenerations
!            IF(iostatus>0) THEN
!                PRINT "(A)", "Write file fail."
!                RETURN
!            END IF
!            WRITE(11, "(F7.3)") stat_scorePerGeneration(i)
!        END DO
!        CLOSE (11, STATUS = 'KEEP')
!
!        OPEN( UNIT=11, FILE='stat_scorePerBestGene.data', ACTION='WRITE', IOSTAT=iostatus )
!        DO i=1, maxGenerations
!            IF(iostatus>0) THEN
!                PRINT "(A)", "Write file fail."
!                RETURN
!            END IF
!            WRITE(11, "(F7.3)") stat_scorePerBestGene(i)
!        END DO
!
!        CLOSE (11, STATUS = 'KEEP')

        DEALLOCATE(stat_scorePerGeneration, stat_scorePerBestGene)
    END SUBROUTINE

    SUBROUTINE stat_mutate(gen, score)
        INTEGER (KIND=integerSize), INTENT(IN)    :: gen
        REAL (KIND=realSize), INTENT(IN)        :: score(popSize)
        REAL (KIND=stat_realSize), SAVE         :: tmp

        !SCORE
        !CROSS
        stat_CrossPerGeneration =stat_CrossPerGeneration + offspringSize

        !MUTATION
        IF(stat_mutationsPerGeneration > stat_maxMutationPerGeneration) stat_maxMutationPerGeneration = stat_mutationsPerGeneration
        IF(stat_mutationsPerGeneration < stat_minMutationPerGeneration) stat_minMutationPerGeneration = stat_mutationsPerGeneration

        stat_averageMutations = stat_averageMutations + REAL(stat_mutationsPerGeneration)
        tmp = stat_averageMutations
        tmp = tmp / REAL(gen)


        PRINT "(A12,I)", "Generation: ", gen

        PRINT "(A)", "Score statistics--------------------------"
        PRINT "(A,F)", "Pop Mean:", mean(score)
        PRINT "(A,F)", "PopKeep Mean:", mean(score(1:popKeep))
        PRINT "(A,F)", "Best gene score:", score(1)

!        PRINT "(A)", "Cross statistics--------------------------"
!        PRINT "(A,I)", "Cross:", offspringSize
!        PRINT "(A,I)", "Total Cross:", stat_CrossPerGeneration
!
!        PRINT "(A)", "Mutations statistics----------------------"
!        PRINT "(A,I)", "Mutations:", stat_mutationsPerGeneration
!        PRINT "(A,F)", "Mean:", tmp
!        PRINT "(A,I)", "Max:", stat_maxMutationPerGeneration
!        PRINT "(A,I)", "Min:", stat_minMutationPerGeneration
!        PRINT *
        PRINT *

        stat_mutationsPerGeneration = 0
    CONTAINS
        REAL (KIND=realSize) FUNCTION mean(data)
            IMPLICIT NONE
            REAL (KIND=realSize), INTENT(IN)    :: data(:)
            INTEGER (KIND=stat_integerSize)     :: i
            mean = 0
            DO i=1, SIZE(data)
                mean = mean + data(i)
            END DO
            mean = mean/SIZE(data)
        END FUNCTION
    END SUBROUTINE


    SUBROUTINE testRandomGeneration
        IMPLICIT NONE !Não permite variaves não declaradas
        INTEGER :: i,j
        INTEGER (KIND=8), PARAMETER :: x = 2000000
        INTEGER (KIND=8), PARAMETER :: y = 500
        REAL :: A(80), B(x)
        REAL :: dummy

        DO i=1 , 80
            A(i) = RAND()
        END DO

        dummy = taskTime(.TRUE.)
        DO i=1 , y
            DO j=1 , x
                B(x) = RAND()
            END DO
        END DO
        dummy = taskTime(.FALSE.)

        PRINT "(A,F7.3,A)","1 billion of random numbers generated in: ", dummy, " [seconds]"

    END SUBROUTINE


    SUBROUTINE printHeader
        IMPLICIT NONE

        NAMELIST /gaparameters/ geneSize, yGeneSize, popSize, maxGenerations, selectionRate,     &
        mutateRate, tournamentSelectionRate, randomSeed, poliCrossSize, resolutionForEncodedReal,&
        popKeep, offspringSize

        WRITE (*, gaparameters)
    END SUBROUTINE

END MODULE



