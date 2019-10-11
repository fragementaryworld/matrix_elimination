program main
    !This program is for getting the solution  of Hx=b with the
    !LU elimination, H is the Hilbert Matrix and b is a
    !column vector whose elements all are one. n is the rank of
    !the Hilbert Matrix. t is the time for solving the n linear
    !equations Hx=b. And the file "lu.dat" will include the time
    !t where n is 10,20,30,...,220.
    implicit none
    integer(kind=4)::n,i
    real(kind=8)::a,b,t
    open(unit=1,file="lu.dat",form="formatted")
    do n = 10,220,10
        call CPU_TIME(a)
        do i=1,100
            call lu(n)
        enddo
        call CPU_TIME(b)
        t = (b-a)/100
        write(unit=1,fmt = 200)n,t
    enddo
    close(1)
    200 format(I4,5X,F30.20)
end program main

subroutine lu(n)
    implicit none
    integer(kind=4)::n,i,j,k
    real(kind=8)::t
    real(kind=8),dimension(n,n)::L,U,H
    real(kind=8),dimension(n)::X,Y,B
    do j=1,n
        do i=1,n
            H(i,j)=1.0D0/(i+j-1)
            U(i,j)=0.0D0
            L(i,j)=0.0D0
            if (i==j) then
                L(i,j)=1.0D0
            endif
        enddo
    enddo
    do i=1,n
        B(i)=1.0D0
    enddo

    do j=1,n
        do i=1,n
            if (i>j) then
                if (j==1) then
                    L(i,j)=H(i,j)/U(j,j)
                else
                    t=0.0D0
                    do k=1,j-1
                        t = t + L(i,k)*U(k,j)
                    enddo
                    L(i,j)=(H(i,j)-t)/U(j,j)
                endif
            else
                if (i==1) then
                    U(i,j)=H(i,j)
                else
                    t=0.0D0
                    do k=1,i-1
                        t = t + L(i,k)*U(k,j)
                    enddo
                    U(i,j)=H(i,j)-t
                endif
            endif
        enddo
    enddo

    do i=1,n
        if (i==1) then
            Y(i)=B(i)
        else
            t=0.0D0
            do k=1,i-1
                t = t + L(i,k)*Y(k)
            enddo
            Y(i)=B(i)-t
        endif
    enddo

    do i=n,1,-1
        if (i==n) then
            X(i)=Y(i)/U(i,i)
        else
            t=0.0D0
            do k=i+1,n
                t = t + U(i,k)*X(k)
            enddo
            X(i)=(Y(i)-t)/U(i,i)
        endif
    enddo

    ! write(*,100)matmul(H,X)
    ! 100 format(4(F30.20,2X))

end subroutine lu