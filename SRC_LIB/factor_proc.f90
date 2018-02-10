SUBROUTINE FACTOR ( N )
    C     Factor N
    C     Test for factors of
    C     2,3,5,... using NEXTPRIME
    C     function to get next
    C     prime.
    INTEGER N, NUMBER, PRIME
    NUMBER = N
    PRIME = 2
    !Start with 2
    10    CONTINUE
    IF
    (NUMBER .EQ. 0) RETURN
    IF
    (
    MOD(NUMBER,PRIME)
    .EQ.
    0
    )
    THEN
    !prime
    !is
    !a
    !factor
    WRITE(6,1)
    PRIME
    NUMBER
    =
    NUMBER
    /
    PRIME
    IF
    (NUMBER
    .EQ.
    1)
    RETURN
ELSE
    PRIME
    =
    NEXTPRIME(PRIME+1)
    !Else
    !test
    !the
    !next
    !prime
ENDIF
GOTO
10
1
FORMAT(1X,'FACTOR
',I)
END

INTEGER
FUNCTION
    NEXTPRIME(
    N
    )
    C     Return
    C     next
    C     prime
    C     number
    C     .GE.
    C     N
    C
    C     For
    C     N
    C     .GE.
    C     1,018,801
    C     returns
    C     next
    C     number
    C     relatively
    C     prime
    C     to
    C     all
    C     primes
    C     below
    C     1009.
    C     Note:
    C     1009
    C     is
    C     next
    C     prime
    C     after
    C     997,
    C     and
    C     1,018,801
    C     =
    C     1009^2.
    C           The
    C           database
    C           may
    C           be
    C           extended
    C           if
    C           larger
    C           primes
    C           are
    C           desired.
    C
    C     Outline:
    C         1.
    C         If
    C         number
    C         .LE.
    C         3
    C         then
    C         return
    C         number.
    C             Number
    C             is
    C             prime,
    C             zero,
    C             or
    C             negative.
    C         2.
    C         If
    C         number
    C         between
    C         4
    C         and
    C         997
    C         then
    C         search
    C         PRIME
    C         array
    C             for
    C             next
    C             prime
    C             .GE.
    C             number.
    C             An
    C             approximating
    C             polynomial
    C             is
    C             used
    C             to
    C             give
    C             us
    C             a
    C             good
    C             initial
    C             guess
    C             for
    C             the
    C             index
    C             I
    C             so
    C             that
    C             PRIME(I)
    C             is
    C             very
    C             close
    C             to
    C             the
    C             prime
    C             we
    C             seek.
    C             The
    C             approximating
    C             polynomial
    C             is:
    C
    C                 I
    C                 =
    C                 -0.3034277E-04*X^2
    C                 +
    C                 0.1918667*X
    C                 +
    C                 8.0918350
    C
    C              The
    C              optimum
    C              coefficients
    C              for
    C              this
    C              polynomial
    C              were
    C              obtained
    C              using
    C              the
    C              OPTIMIZE
    C              program
    C              in
    C              the
    C              [.OPTIMIZE]
    C              directory.
    C              Look
    C              there
    C              for
    C              further
    C              details.
    C              Thus,
    C              this
    C              function
    C              is
    C              right
    C              on
    C              the
    C              money
    C              452
    C              times,
    C              is
    C              low
    C              by
    C              one
    C              224
    C              times,
    C              is
    C              low
    C              by
    C              two
    C              40
    C              times,
    C              and
    C              is
    C              low
    C              by
    C              three
    C              1
    C              time.
    C              It's
    C              high
    C              by
    C              one
    C              185
    C              times,
    C              high
    C              by
    C              two
    C              68
    C              times,
    C              high
    C              by
    C              three
    C              12
    C              times,
    C              high
    C              by
    C              four
    C              7
    C              times,
    C              high
    C              by
    C              five
    C              4
    C              times,
    C              and
    C              high
    C              by
    C              6
    C              one
    C              time.
    C
    C         3.
    C         Else
    C         number
    C         is
    C         greater
    C         than
    C         997.
    C         Do
    C         normal
    C         search
    C         for
    C         next
    C         prime.
    C             3a.
    C             Increment
    C             NUMBER
    C             to
    C             an
    C             odd
    C             number.
    C             Since
    C             even
    C             numbers
    C             are
    C                  not
    C                  prime
    C                  we
    C                  skip
    C                  over
    C                  even
    C                  numbers.
    C             3b.
    C             Test
    C             if
    C             NUMBER
    C             is
    C             evenly
    C             divisible
    C             by
    C             any
    C             prime
    C             in
    C             the
    C                  PRIME
    C                  array
    C                  less
    C                  than
    C                  SQRT(NUMBER).
    C                  If
    C                  NUMBER
    C                  MOD
    C                  PRIME(I)
    C                  =
    C                  0,
    C                  then
    C                  NUMBER
    C                  is
    C                  not
    C                  prime.

    INTEGER
    N,
    NUMBER
    INTEGER
    NUMPRIMES
    PARAMETER
    (NUMPRIMES
    =
    168)
    INTEGER
    PRIME(NUMPRIMES)
    DATA
    PRIME
    /
    *   2,
    *   3,
    *   5,
    *   7,
    *   11,
    *   13,
    *   17,
    *   19,
    *   23,
    *   29,
    *   31,
    *   37,
    *   41,
    *  43,
    *  47,
    *  53,
    *  59,
    *  61,
    *  67,
    *  71,
    *  73,
    *  79,
    *  83,
    *  89,
    *  97,
    *  101,
    * 103,
    * 107,
    * 109,
    * 113,
    * 127,
    * 131,
    * 137,
    * 139,
    * 149,
    * 151,
    * 157,
    * 163,
    * 167,
    * 173,
    * 179,
    * 181,
    * 191,
    * 193,
    * 197,
    * 199,
    * 211,
    * 223,
    * 227,
    * 229,
    * 233,
    * 239,
    * 241,
    * 251,
    * 257,
    * 263,
    * 269,
    * 271,
    * 277,
    * 281,
    * 283,
    * 293,
    * 307,
    * 311,
    * 313,
    * 317,
    * 331,
    * 337,
    * 347,
    * 349,
    * 353,
    * 359,
    * 367,
    * 373,
    * 379,
    * 383,
    * 389,
    * 397,
    * 401,
    * 409,
    * 419,
    * 421,
    * 431,
    * 433,
    * 439,
    * 443,
    * 449,
    * 457,
    * 461,
    * 463,
    * 467,
    * 479,
    * 487,
    * 491,
    * 499,
    * 503,
    * 509,
    * 521,
    * 523,
    * 541,
    * 547,
    * 557,
    * 563,
    * 569,
    * 571,
    * 577,
    * 587,
    * 593,
    * 599,
    * 601,
    * 607,
    * 613,
    * 617,
    * 619,
    * 631,
    * 641,
    * 643,
    * 647,
    * 653,
    * 659,
    * 661,
    * 673,
    * 677,
    * 683,
    * 691,
    * 701,
    * 709,
    * 719,
    * 727,
    * 733,
    * 739,
    * 743,
    * 751,
    * 757,
    * 761,
    * 769,
    * 773,
    * 787,
    * 797,
    * 809,
    * 811,
    * 821,
    * 823,
    * 827,
    * 829,
    * 839,
    * 853,
    * 857,
    * 859,
    * 863,
    * 877,
    * 881,
    * 883,
    * 887,
    * 907,
    * 911,
    * 919,
    * 929,
    * 937,
    * 941,
    * 947,
    * 953,
    * 967,
    * 971,
    * 977,
    * 983,
    * 991,
    * 997
    * /

    NUMBER
    =
    N

    IF
    (NUMBER
    .LE.
    3)
    THEN
    !if
    !number
    !.le.
    !3
    !then
    !return
    !number
    NEXTPRIME
    =
    NUMBER
    RETURN
ENDIF

IF
(NUMBER
.LE.
997)
THEN
!just
!search
!thru
!PRIME
!array
I
=
(-0.3034277E-04*NUMBER
+
0.1918667)*NUMBER
+
8.0918350
!(see
!comments
!above)
I
=
MIN(I,NUMPRIMES)
!don't
!let
!I
!go
!above
!997
DO
WHILE
(PRIME(I)
.LT.
NUMBER)
!Search
!upward
!for
!first
!prime
!greater
!than
!NUMBER
I
=
I
+
1
ENDDO
DO
WHILE
(PRIME(I)
.GE.
NUMBER)
!Search
!downward
!for
!first
!prime
!less
!than
!NUMBER,
I
=
I
-
1
ENDDO
NEXTPRIME
=
PRIME(I+1)
!then
!take
!next
!prime.
RETURN

  ELSE
      !Else
      !normal
      !search
      !for
      !next
      !prime

      IF
      (
      MOD(NUMBER,2)
      .EQ.
      0
      )
      NUMBER
      =
      NUMBER
      +
      1
      !rule
      !out
      !even
      !numbers.
      !They
      !are
      !not
      !prime
      I
      =
      2
      !start
      !with
      !second
      !prime
      !since
      !we
      !only
      !test
      !odd
      !numbers.
      DO
      WHILE
      (
      PRIME(I)*PRIME(I)
      .LE.
      NUMBER
      )
      !test
      !for
      !all
      !primes
      !.LE.
      !SQRT(NUMBER)
      IF
      (
      MOD(NUMBER,PRIME(I))
      .EQ.
      0
      )
      THEN
      !it's
      !not
      !prime
      NUMBER
      =
      NUMBER
      +
      2
      !even
      !numbers
      !aren't
      !prime
      I
      =
      2
      !start
      !over
      !again
      !with
      !new
      !number
  ELSE
      I
      =
      I
      +
      1
      !test
      !with
      !next
      !prime
      IF
      (
      I
      .GT.
      NUMPRIMES)
      GOTO
      10
      !exit
      !loop.
      !NUMBER
      !is
      !relatively
      !prime
      !to
      !first
      !168
      !primes.
  ENDIF
  ENDDO
  10
  NEXTPRIME
  =
  NUMBER
  RETURN

  ENDIF

  END
  !! Local Variables:
  !! mode: f90
  !! show-trailing-whitespace: t
  !! coding: utf-8
  !! f90-do-indent: 4
  !! f90-if-indent: 4
  !! f90-type-indent: 4
  !! f90-program-indent: 4
  !! f90-continuation-indent: 4
  !! End:
  !! vim: set sw=4 ts=4 et tw=80 smartindent :

