function dist(a, b)
use constant
real(dp) :: dist
real(dp), intent(in) :: a(3), b(3)
real(dp) :: tmp(3)

    tmp = a - b
    dist = dot_product( tmp, tmp )
    dist = dsqrt( dist )

end function

! subroutine : law of cosine
! input - lenght of each side a, b, c
! output - the angle of C gamma in rad
subroutine cosineAngle(a,b,c, gamma)
use constant
implicit none
real(dp), intent(in) :: a,b,c
real(dp), intent(out) :: gamma

    gamma = (a*a + b*b - c*c) / (2d0*a*b)
    if ( gamma .gt. 1d0 ) then
        write(*,*) 'cos(alpha) is greater than 1! '
        write(*,*) a,b,c
        stop
    end if
    
    gamma = dacos(gamma)
    !write(*,*) 'gamma: ', gamma / pi * 180d0

end subroutine