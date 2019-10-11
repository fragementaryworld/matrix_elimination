program main
    !This program is for getting the solution  of Hx=b with the
    !Gauss-Jordan elimination, H is the Hilbert Matrix and b is a
    !column vector whose elements all are one. n is the rank of
    !the Hilbert Matrix. t is the time for solving the n linear
    !equations Hx=b. And the file "gj.dat" will include the time
    !t where n is 200,400,600,...,3000.

    implicit none
    integer(kind=4)::n
    real(kind=8)::a,b,t
    open(unit=1,file="gj.dat",form="formatted")
    do n = 200,3000,200
        call CPU_TIME(a)
        call gj(n)
        call CPU_TIME(b)
        t = (b-a)
        write(unit=1,fmt = 200)n,t
    enddo
    close(1)
    200 format(I4,5X,F30.20)
end program main

subroutine gj(n)
    implicit none
    integer(kind=4)::n,i,j,k
    real(kind=8)::d,t
    real(kind=8),dimension(n,n)::H
    real(kind=8),dimension(n,n+1)::C
    real(kind=8),dimension(n)::X
    integer(kind=4),dimension(n)::js,is

    do j = 1,n
        do i=1,n
            C(i,j)=1.0D0/(i+j-1)
        enddo
    enddo
    do i=1,n
        C(i,n+1)=1.0D0
    enddo

    do k=1,n
        d = 0.0D0
        do j=k,n
            do i=k,n
                if (abs(C(i,j)) > d) then
                    d = abs(C(i,j))
                    is(k) = i
                    js(k) = j
                endif
            enddo
        enddo
        if (abs(d)==0.0D0) then
            write(*,*)"The matrix is singular!"
        endif
        if (is(k) /= k) then
            do j=1,n+1
                t = C(is(k),j)
                C(is(k),j) = C(k,j)
                C(k,j) = t
            enddo
        endif
        if (js(k) /= k) then
            do i=1,n
                t = C(i,js(k))
                C(i,js(k)) = C(i,k)
                C(i,k) = t
            enddo
        endif
        do j=1,n+1
            C(k,j) = C(k,j)/d
        enddo
        do i=1,n
            if (i /= k) then
                t = C(i,k)
                do j=1,n+1
                    C(i,j) = C(i,j) - C(k,j)*t
                enddo
            endif
        enddo
    enddo
    X = C(:,n+1)
    do k=n,1,-1
        t = X(k)
        X(k) = X(js(k))
        X(js(k)) = t
    enddo

    ! do j=1,n
    !     do i=1,n
    !         H(i,j) = 1.0D0/(i+j-1)
    !     enddo
    ! enddo
    ! write(*,100)matmul(H,X)
    ! 100 format(4(F30.20,2X))
end subroutine gj