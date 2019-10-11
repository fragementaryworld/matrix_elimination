program main
    !This program is for getting the solution  of Hx=b with the
    !Gauss elimination, H is the Hilbert Matrix and b is a
    !column vector whose elements all are one. n is the rank of
    !the Hilbert Matrix. t is the time for solving the n linear
    !equations Hx=b. And the file "lu.dat" will include the time
    !t where n is 10,20,30,...,220.
    implicit none
    integer(kind=4)::n,i
    real(kind=8)::a,b,t
    open(unit=1,file="gauss.dat",form="formatted")
    do n = 10,220,10
        call CPU_TIME(a)
        do i =1,100
            call guass(n)
        enddo
        call CPU_TIME(b)
        t = (b-a)/100
        write(unit=1,fmt = 200)n,t
    enddo
    close(1)
    200 format(I4,5X,F30.20)
end program main

subroutine guass(n)
    implicit none
    integer(kind=4)::n
    real(kind=8),dimension(n,n+1)::C
    real(kind=8),dimension(n,n)::H
    real(kind=8),dimension(n)::x
    real(kind=8) :: t
    integer(kind=4)::i,j,k

    do j = 1,n
        do i=1,n
            C(i,j)=1.0D0/(i+j-1)
        enddo
    enddo
    do i=1,n
        C(i,n+1)=1.0D0
    enddo

    do k=1,n
        t = C(k,k)
        do j=1,n+1
            C(k,j)=C(k,j)/t
        enddo
        do i=1,n
            if (i /= k) then
                t = C(i,k)
                do j=1,n+1
                    C(i,j)=C(i,j)-C(k,j)*t
                enddo
            endif
        enddo
    enddo

    x=C(:,n+1)
    ! write(*,100)x

    ! do j = 1,n
    !     do i = 1,n
    !         H(i,j)=1.D0/(i+j-1)
    !     enddo
    ! enddo
    ! write(*,100)matmul(H,x)
    ! 100 format(4(F30.20,2X))
end subroutine