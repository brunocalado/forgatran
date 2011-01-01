MODULE parameters
    IMPLICIT NONE
!CONSTANTS----------------------------------------------------------------------
    INTEGER (KIND=1), PARAMETER :: logicalSize     = 1
    INTEGER (KIND=1), PARAMETER :: characterSize   = 1
    INTEGER (KIND=1), PARAMETER :: integerSize     = 8
    INTEGER (KIND=1), PARAMETER :: integerGeneSize = 1
    INTEGER (KIND=1), PARAMETER :: realSize        = 8
    INTEGER (KIND=1), PARAMETER :: continousSize   = 4

!DEBUG--------------------------------------------------------------------------
!DEC$ DEFINE debug = 0

!-------------------------------------------------------------------------------

    INTEGER (KIND=integerSize) :: geneSize      = 100    !> 1D gene size. Recalculated
    INTEGER (KIND=integerSize) :: yGeneSize     = 1     !> 2D gene size. Number of parameters for real or integer encode
    INTEGER (KIND=integerSize) :: totalGeneSize                   !>Self-calculated
    INTEGER (KIND=integerSize) :: popSize       = 20    !> Population size. Lower pop allowed 2

    !>Total of generations = maxGenerations x extendsGenerations
    INTEGER (KIND=integerSize) :: maxGenerations = 10000

    REAL (KIND=realSize) :: selectionRate = 0.500D0     !> Piece of population selected for mate
    REAL (KIND=realSize) :: offspringRate = 0.500D0     !> Amount of offspring produced

    REAL (KIND=realSize) :: mutateRate           = 0.0500D0    !> Starting mutate probability per bit
    REAL (KIND=realSize) :: highestMutateRate    = 0.0500D0    !> Maximum mutate probability per bit
    REAL (KIND=realSize) :: lowestMutateRate     = 0.0010D0    !> Minimum mutate probability per bit
    REAL (KIND=realSize) :: mutateRateMultiplier = 1.2500D0    !>
    REAL (KIND=realSize) :: highestDelta         = 0.1100D0    !>
    REAL (KIND=realSize) :: lowestDelta          = 0.0500D0    !>

    LOGICAL (KIND=logicalSize) :: insertion               = .FALSE.
    REAL (KIND=realSize)       :: tournamentSelectionRate = 0.7000D0    !> Probability to selected best gene

    INTEGER (KIND=integerSize) :: randomSeed    = 1003      !>Set to -1 for random choice using date and time
    INTEGER (KIND=integerSize) :: poliCrossSize = 10         !>LOWER SIZE 3

    !> INTEGER Array encode parameters
!    INTEGER (KIND=integerSize), ALLOCATABLE        :: upperBoundForEncodedInteger(:)
!    INTEGER (KIND=integerSize), ALLOCATABLE        :: lowerBoundForEncodedInteger(:)

    INTEGER (KIND=integerSize) :: popKeep                              !>Self-calculated
    INTEGER (KIND=integerSize) :: offspringSize                        !>Self-calculated

!EXCLUSIVE FOR BINARY ENCODE----------------------------------------------------
    !> REAL Array encode parameters
    INTEGER (KIND=integerSize)        :: resolutionForEncodedReal    = 1024*16     !> Must be 2**n, where n is a integer
    REAL (KIND=realSize), ALLOCATABLE :: upperBoundForEncodedReal(:)
    REAL (KIND=realSize), ALLOCATABLE :: lowerBoundForEncodedReal(:)
    REAL (KIND=realSize), ALLOCATABLE :: incrementForEncodedReal(:)                !>Self-calculated


!EXCLUSIVE FOR INTEGER ENCODE---------------------------------------------------
    INTEGER (KIND=integerSize) :: minimumInteger = 1
    INTEGER (KIND=integerSize) :: maximumInteger = 10


!EXCLUSIVE FOR REAL ENCODE------------------------------------------------------


!POPULATION---------------------------------------------------------------------
    LOGICAL (KIND=logicalSize), ALLOCATABLE :: pop(:, :)
    LOGICAL (KIND=logicalSize), ALLOCATABLE :: popAux(:, :)
!
!    INTEGER (KIND=integerGeneSize), ALLOCATABLE :: pop(:, :)
!    INTEGER (KIND=integerGeneSize), ALLOCATABLE :: popAux(:, :)

!    REAL (KIND=continousSize), ALLOCATABLE :: pop(:, :)
!    REAL (KIND=continousSize), ALLOCATABLE :: popAux(:, :)


!DON'T TOUCH THESE--------------------------------------------------------------
    LOGICAL (KIND=logicalSize), ALLOCATABLE :: mask(:)
    LOGICAL (KIND=logicalSize), ALLOCATABLE :: maskNiche(:)

    INTEGER (KIND=integerSize), ALLOCATABLE :: geneForMate(:)
    REAL (KIND=realSize), ALLOCATABLE       :: rankweighProbability(:)
    REAL (KIND=realSize), ALLOCATABLE       :: costWeighProbability(:)
    REAL (KIND=realSize), ALLOCATABLE       :: geneScore(:)



END MODULE
