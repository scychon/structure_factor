!****************************************************************
! 
! this module contains utility functions for vector calculation
!
!****************************************************************

MODUlE fvector
    implicit none
 
contains

    ! Subroutines to expand the size of allocatable arrays
    subroutine expand1D ( array, nsize )
        implicit none
        integer, intent(in) :: nsize
        real*8, intent(inout), allocatable, dimension(:) :: array
        real*8, allocatable, dimension(:) :: temp
        integer :: isize
        
        isize = size(array(:))
        allocate(temp(isize+nsize))
        temp=0
        temp(:isize)=array
        deallocate(array)
        allocate(array(isize+nsize))
        array=temp
        deallocate(temp)
    end subroutine expand1D

    subroutine expand2D ( array, nsize1, nsize2 )
        implicit none
        integer, intent(in) :: nsize1, nsize2
        real*8, intent(inout), allocatable, dimension(:,:) :: array
        real*8, allocatable, dimension(:,:) :: temp
        integer :: isize1, isize2,i
        
        isize1 = size(array(:,1))
        isize2 = size(array(1,:))
        allocate(temp(isize1+nsize1,isize2+nsize2))
        temp=0
        temp(:isize1,:isize2)=array
        deallocate(array)
        allocate(array(isize1+nsize1,isize2+nsize2))
        array=temp
        deallocate(temp)
    end subroutine expand2D

    subroutine expand3D ( array, nsize1, nsize2, nsize3 )
        implicit none
        integer, intent(in) :: nsize1, nsize2, nsize3
        real*8, intent(inout), allocatable, dimension(:,:,:) :: array
        real*8, allocatable, dimension(:,:,:) :: temp
        integer :: isize1, isize2, isize3,i
        
        isize1 = size(array(:,1,1))
        isize2 = size(array(1,:,1))
        isize3 = size(array(1,1,:))
        allocate(temp(isize1+nsize1,isize2+nsize2,isize3+nsize3))
        temp=0

        temp(:isize1,:isize2,:isize3)=array
        deallocate(array)
        allocate(array(isize1+nsize1,isize2+nsize2,isize3+nsize3))
        array=temp
        deallocate(temp)
    end subroutine expand3D


    ! Subroutines to schrink the size of allocatable arrays to fit the dimension
    subroutine shrink1D ( array, nsize )
        implicit none
        integer, intent(in) :: nsize
        real*8, intent(inout), allocatable, dimension(:) :: array
        real*8, allocatable, dimension(:) :: temp
        
        allocate(temp(nsize))
        temp=0
        temp(:)=array(:nsize)
        deallocate(array)
        allocate(array(nsize))
        array=temp
        deallocate(temp)
    end subroutine shrink1D

    subroutine shrink2D ( array, nsize1, nsize2 )
        implicit none
        integer, intent(in) :: nsize1, nsize2
        real*8, intent(inout), allocatable, dimension(:,:) :: array
        real*8, allocatable, dimension(:,:) :: temp
        
        allocate(temp(nsize1,nsize2))
        temp=array(:nsize1,:nsize2)
        deallocate(array)
        allocate(array(nsize1,nsize2))
        array=temp
        deallocate(temp)
    end subroutine shrink2D

    subroutine shrink3D ( array, nsize1, nsize2, nsize3 )
        implicit none
        integer, intent(in) :: nsize1, nsize2, nsize3
        real*8, intent(inout), allocatable, dimension(:,:,:) :: array
        real*8, allocatable, dimension(:,:,:) :: temp
        
        allocate(temp(nsize1,nsize2,nsize3))
        temp=array(:nsize1,:nsize2,:nsize3)
        deallocate(array)
        allocate(array(nsize1,nsize2,nsize3))
        array=temp
        deallocate(temp)
    end subroutine shrink3D

    function vecAdd( vec1, vec2 )
        implicit none
        real*8, intent(in), dimension(3) :: vec1, vec2
        real*8, dimension(3) :: vecAdd
        vecAdd(:) = vec1(:) - vec2(:)
    end function vecAdd

    function vecCross( vec1, vec2 )
        implicit none
        real*8, intent(in), dimension(3) :: vec1, vec2
        real*8, dimension(3) :: vecCross
        vecCross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        vecCross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        vecCross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    end function vecCross

    function vecUnitNorm( vec1, vec2 )
        implicit none
        real*8, intent(in), dimension(3) :: vec1, vec2
        real*8, dimension(3) :: vecUnitNorm
        real*8, dimension(3) :: vecNorm
        vecNorm = vecCross(vec1,vec2)
        vecUnitNorm(:) = vecNorm(:)/norm(vecNorm)
    end function vecUnitNorm
    

    real*8 function getdr( vec_dr, vec_box )
        implicit none
        real*8, intent(in), dimension(3) :: vec_box, vec_dr
        real*8, dimension(3) :: tempvec
    
        tempvec(:) = vec_dr(:) - (nint(vec_dr(:)/vec_box(:)))*vec_box(:)
        getdr = norm(tempvec)
    end function getdr
    
    real*8 function norm( vec )
        implicit none
        real*8, intent(in), dimension(3) :: vec
    
        norm = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
    
    end function norm
 
end module fvector
