program main
integer i

i = 1
call random_seed()
do while (.true.)
write(*,*) i
call testBonCar
i = i+1
end do

end program


subroutine testBonCar
use constant
implicit none
integer :: Natoms, Nbonds
real(dp), allocatable :: q(:,:), b(:), q1(:,:), b1(:)
integer :: res ! result
integer :: i

    open(999, file='failure.xyz', status='old', action='write')
    Natoms = 10
    Nbonds = (Natoms - 1) * Natoms / 2
    
    allocate( q(3,Natoms), b(Nbonds) )
    allocate( q1(3,Natoms), b1(Nbonds) )
    
    call genRandCar(Natoms, q)
    call calcBond(Natoms, q, Nbonds, b)
    call bonCar(Nbonds, b, Natoms, q1, res)
    call calcBond(Natoms, q1, Nbonds, b1)
    !write(*,*) 'Cart:'
    !write(*,'(3F12.6)') q1
    !write(*,*) 'bond list: '
    !write(*,'(6F12.6)') b1
    !write(*,*) 'Cartesian difference : '
    !write(*,'(3E12.2)') q1 - q
    !write(*,*) 'bond lengths difference : '
    if ( res .lt. 0 ) then
        write(999,*) Natoms
        write(999,*) 'Wenbin, FAN'
        do i = 1, Natoms
            write(999,'(A2,3E25.15)') ' H', q(1:3,i)
        end do
        write(*,*) 'original Cartesian: '
        write(*,'(3F12.6)') q
        write(*,*) 'original bond list: '
        write(*,'(6F12.6)') b
        write(*,*) 'Cart:'
        write(*,'(3F12.6)') q1
        write(*,*) 'bond list: '
        write(*,'(6F12.6)') b1
        write(*,*) 'Cartesian difference : '
        stop
    end if

end subroutine


subroutine bonCar(Nbonds, b, N, q, res)
use constant
implicit none
integer, intent(in) :: Nbonds, N
real(dp), intent(in) :: b(Nbonds)
real(dp), intent(out) :: q(3, N)
integer :: res ! result

integer :: i,j,k, order
real(dp) :: alpha ! angle
real(dp) :: q4(3,2) ! two cart coordinates of 4th point
real(dp) :: bc(3) ! bonds to calculate the 4th point
!real(dp) :: bq(N-4) ! check the orientation of 4th point, `q` means check (cheque -> que -> q)
real(dp), external :: dist
real(dp) :: checkDist
integer :: checkPass
real(dp) :: bout(Nbonds)

    if (Nbonds .le. 0) stop 'Wrong number of bonds! '
    if (N .ne. int( (1d0 + dsqrt(1d0+8d0*Nbonds)) / 2d0 ) ) then
        stop 'Wrong number of atoms! '
    end if
    
    ! check colinear condition
    call colinearCheck(Nbonds, b, N)
    
    ! initial the Cartesian coordinates, set to a huge negetive value
    if (N .le. 0) then
        stop 'Negative number of atoms! '
    else
        q = -huge(1d0)
    end if
    
    ! the first atom was placed at the origin
    if (N .ge. 1) then
        q(1:3,1) = 0d0
    else
        return
    end if
    
    ! the second atom was placed at the x-axis
    if (N .ge. 2) then
        q(1,2) = b(1)
        q(2:3,2) = 0d0
    else
        return
    end if
    
    ! the third atom was placed at the xOy plane
    if (N .ge. 3) then
        call cosineAngle(b(1),b(2),b(N), alpha)
        q(1,3) = b(2)*dcos(alpha)
        q(2,3) = b(2)*dsin(alpha)
        q(3,3) = 0d0
    else
        return
    end if
    
    ! the fourth atom
    !i = 4
    !if (N .eq. 4) then
    !    bc(1) = b(i-1)           !1-4, 3rd
    !    bc(2) = b(N-1 + 2)       !2-4
    !    bc(3) = b(N-1 + N-2 + 1) !3-4
    !    !write(*,*) N-1, N-1 + i-2, N-1 + N-2 + i-3
    !    call calcFourthCart(q(1:3,1:3), bc(1:3), q4(1:3,1:2))
    !    q(1:3,4) = q4(1:3,1)
    !else 
    !    return
    !end if
    
    if (N .ge. 4) then
    do i = 4, N
    
        bc(1) = b(i-1)
        bc(2) = b(N-1 + i-2)
        bc(3) = b(N-1 + N-2 + i-3)

        call calcFourthCart(q(1:3,1:3), bc(1:3), q4(1:3,1:2))
        if(i.eq.4) q(1:3,i) = q4(1:3,1)
        
        ! check and choose the direction of the fourth coordinate
        do j = 1, i-4
            order = ( - (j+2)*(j+2) + 2*(j+2)*N - (j+2) ) / 2 + i-3 - j
            !       \sum_{i=1}^{2+j} (N-i)                    +(i-3)- j
            
            checkPass = -100
            do k = 1, 2
                checkDist = abs( dist(q(1:3,j+3),q4(1:3,k)) - b(order) )
                
                if ( checkDist .lt. eps ) then
                    q(1:3,i) = q4(1:3,k)
                    checkPass = 100
                    exit
                end if
            end do
            if (checkPass .lt. 0) then
                !write(*,*) i, q(1:3,i)
                res = -100
                return
            end if
        end do
        
    end do
    end if
    res = 100
    
    call calcBond(N, q, Nbonds, bout)
    if ( maxval(abs( bout - b )) .gt. 1d-6 ) res = -100
    

end subroutine


subroutine colinearCheck(Nbonds, b, Natoms)
use constant
implicit none
integer, intent(in) :: Nbonds, Natoms
real(dp), intent(in) :: b(Nbonds)

integer i,j,k
integer, external :: getOrder
integer x,y,z, res

    do i = 1, Natoms-2
    do j = i+1, Natoms-1
    do k = j+1, Natoms
        x = getOrder(Natoms, i, j)
        y = getOrder(Natoms, i, k)
        z = getOrder(Natoms, j, k)
        call deltaCheck(b(x), b(y), b(z), res)
        if (res .lt. 0) then
            write(*,'(3I3,A)') i,j,k,' colinear'
        end if
    end do
    end do
    end do

end subroutine


function getOrder(N, a, b)
integer :: getOrder
integer, intent(in) :: N, a, b

    if ( a .le. 0 ) stop 'wrong input number! '
    getOrder = (a-1) * (2*N-a) / 2 + b-a

end function


subroutine deltaCheck(a,b,c,res)
use constant
implicit none
real(dp), intent(in) :: a,b,c
integer, intent(out) :: res

real(dp) :: g(3)

    res = 100
    !if ( a + b .le. c ) then
    !    res = -100
    !end if
    !if ( a + c .le. b ) then
    !    res = -100
    !end if
    !if ( b + c .le. a ) then
    !    res = -100
    !end if
    call cosineAngle(a,b,c, g(1))
    call cosineAngle(b,c,a, g(2))
    call cosineAngle(c,a,b, g(3))
    if (minval(g) .lt. pi/1800d0 ) res = -100

end subroutine


! 2020-12-02 13:53:42 Wenbin, FAN @ SHU
! subroutine : calculate the Cartesian coordinate of fourth point
! source : https://math.stackexchange.com/questions/2969363/finding-a-4th-point-in-3d-space-knowing-3-other-points-and-2-distances-to-the-4t
! input
!   q - three known Cartesian
!   b - three bond lengths from unknown point to three known point
! output
!   q4 - two Cartesian coordinates of the fourth point, one above and one below 123 plane
subroutine calcFourthCart(q, b, q4)
use constant
implicit none
real(dp), intent(in) :: q(3,3), b(3)
real(dp), intent(out) :: q4(3,2)

integer :: i
real(dp) :: d1, d2, d3
real(dp) :: x1,y1,z1, x2,y2,z2, x3,y3,z3
real(dp) :: ex1,ey1,ez1, ex2,ey2,ez2, ex3,ey3,ez3
real(dp) :: hh, ii, jj, tt
real(dp) :: u,v,w, x,y,z

    d1 = b(1); d2 = b(2); d3 = b(3)
    x1 = q(1,1); y1 = q(2,1); z1 = q(3,1)
    x2 = q(1,2); y2 = q(2,2); z2 = q(3,2)
    x3 = q(1,3); y3 = q(2,3); z3 = q(3,3)
    
    ex1 = x2 - x1; ey1 = y2 - y1; ez1 = z2 - z1
    hh = dsqrt( ex1*ex1 + ey1*ey1 + ez1*ez1 )
    if (hh .lt. eps) stop 'too close'
    ex1 = ex1/hh; ey1 = ey1/hh; ez1 = ez1/hh
    !write(*,*) ex1, ey1, ez1
    
    ii = ex1 * (x3 - x1) + ey1 * (y3 - y1) + ez1 * (z3 - z1)
    
    ex2 = x3 - x1 - ii * ex1
    ey2 = y3 - y1 - ii * ey1
    ez2 = z3 - z1 - ii * ez1
    tt = dsqrt( ex2*ex2 + ey2*ey2 + ez2*ez2 )
    if (tt .lt. eps) stop 'colinear'
    ex2 = ex2/tt; ey2 = ey2/tt; ez2 = ez2/tt
    
    jj = ex2 * ( x3 - x1 ) + ey2 * (y3 - y1) + ez2 * (z3 - z1)
    if ( abs(jj) .lt. eps ) stop 'colinear second'
    
    ex3 = ey1*ez2 - ez1*ey2
    ey3 = ez1*ex2 - ex1*ez2
    ez3 = ex1*ey2 - ey1*ex2
    
    u = (d1*d1 - d2*d2 + hh*hh) / (2d0*hh)
    v = (d1*d1 - d3*d3 + ii*(ii - 2d0*u) + jj*jj) / (2d0*jj)
    w = d1*d1 - u*u - v*v
    
    if ( w .lt. 0d0 ) stop 'no solution for fourth point'
    w = dsqrt(w)
    !if ( w .lt. eps ) then
    !    write(*,*) 'small w', w
    !end if
    !    x = x1 + u * ex1 + v * ex2
    !    y = y1 + u * ey1 + v * ey2
    !    z = z1 + u * ez1 + v * ez2
    !else
    !    x = x1 + u * ex1 + v * ex2 + w * ex3
    !    y = y1 + u * ey1 + v * ey2 + w * ey3
    !    z = z1 + u * ez1 + v * ez2 + w * ez3
    !end if
    q4(1,1) = x1 + u * ex1 + v * ex2 + w * ex3
    q4(2,1) = y1 + u * ey1 + v * ey2 + w * ey3
    q4(3,1) = z1 + u * ez1 + v * ez2 + w * ez3
    q4(1,2) = x1 + u * ex1 + v * ex2 - w * ex3
    q4(2,2) = y1 + u * ey1 + v * ey2 - w * ey3
    q4(3,2) = z1 + u * ez1 + v * ez2 - w * ez3

end subroutine


subroutine calcBond(Natoms, q, Nbonds, b)
use constant
implicit none
integer, intent(in) :: Natoms, Nbonds
real(dp), intent(in) :: q(3, Natoms)
real(dp), intent(out) :: b(Nbonds)
integer :: i,j,k
real(dp) :: tmp(3) ! temporary vector

    if ( Natoms .le. 1 ) stop 'Wrong number of atoms! '
    if ( Nbonds .ne. (Natoms - 1) * Natoms / 2 ) stop 'Wrong number of bonds! '
    
    k = 0
    do i = 1, Natoms - 1
        do j = i+1, Natoms
            k = k + 1
            tmp(1:3) = q(1:3,j) - q(1:3,i)
            b(k) = dot_product( tmp(1:3), tmp(1:3) )
            b(k) = dsqrt( b(k) )
        end do
    end do
    if ( k .ne. Nbonds ) stop 'Wrong bonds results! '

end subroutine


subroutine genRandCar(Natoms, q)
use constant
implicit none
integer, intent(in) :: Natoms
real(dp), intent(out) :: q(3, Natoms)
real(dp), parameter :: s = 5d0 ! scale

    call random_number( q )
    q = 0d0
    if (Natoms .eq. 0) stop 'The number of atom can not be zero! '
    if (Natoms .ge. 1) q(:,1) = 0d0
    
    if (Natoms .ge. 2) then
        call random_number( q(1,2) ) ! atom 2 is placed on the x axis
        q(1,2) = ( q(1,2) - 0.5d0 ) * s
    end if
    
    if (Natoms .ge. 3) then
        call random_number( q(1:2,3) ) ! atom 3 is placed on the xOy plane
        q(1:2,3) = ( q(1:2,3) - 0.5d0 ) * s
    end if
    
    if (Natoms .ge. 4) then
        call random_number( q(:,4:Natoms) )
        q(:,4:Natoms) = q(:,4:Natoms) - 0.5d0
        q(:,4:Natoms) = q(:,4:Natoms) * s
    end if

end subroutine