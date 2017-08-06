program testfile
    implicit none
    integer :: i
    integer, dimension(2,2,2) :: box
    integer, dimension(2,2) :: box3
    integer, dimension(3) :: box2
    box = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8 /), (/2,2,2/))
    !box2 = (/ (box(i,i), i = 1,3) /)

    write(*,*) box(:,1,1)
    write(*,*) box(:,1,2)
    write(*,*) box
    write(*,*) sum(box,dim=1)
    box3 = sum(box,dim=1)
    write(*,*) box3(:,1)
    write(*,*) box3(:,2)
    write(*,*) 'sum(box3(:,1:2),dim=1)'
    write(*,*) sum(box3(:,1:1),dim=1)
end program
