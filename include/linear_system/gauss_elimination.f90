!===============================================================
! Module: gauss_elimination
! Purpose: Solve linear system Ax = b using guass decomposition
!          with pving
!===============================================================

module gauss_elimination
  implicit none

contains
  subroutine gauss_eliminate(A, B, x, N, err)
    integer, intent(in)     :: N
    real(8), intent(inout)  :: A(N, N)
    real(8), intent(inout)  :: B(N)
    real(8), intent(out)    :: x(N)
    integer, intent(out)    :: err
    real(8), parameter :: TOL = 1.0e-10
    real(8) :: coef, temp, max_val
    integer :: pv
    integer :: i, j

    ! Forward elimination
    do i = 1, N
      max_val = abs(A(i,i))
      pv = i
      do j = i+1, N
        if (abs(A(j,i)) > max_val) then
          max_val = abs(A(j,i))
          pv = j
        end if
      end do

      if (max_val < TOL) then
        err = 1
        return
      end if

      if (pv /= i) then
        A([i,pv],:) = A([pv,i],:)
        B([i,pv]) = B([pv,i])
      end if

      do j = i+1, N
        coef = A(j,i) / A(i,i)
        A(j, i+1:N) = A(j, i+1:N) - coef * A(i, i+1:N)
        B(j) = B(j) - coef * B(i)
        A(j, i) = 0.0d0
      end do
    end do

    ! Back substitution
    do i = N, 1, -1
      x(i) = (B(i) - sum(A(i,i+1:N)*x(i+1:N))) / A(i,i)
    end do
  end subroutine gauss_eliminate
    
  subroutine gauss_eliminate_fpv(A, B, x, N, err)
    integer, intent(in)     :: N
    real(8), intent(inout)  :: A(N, N)
    real(8), intent(inout)  :: B(N)
    real(8), intent(out)    :: x(N)
    integer, intent(out)    :: err
    real(8), parameter :: TOL = 1.0d-10
    real(8) :: coef, temp_vect(N), temp, max_val
    integer :: col_pv(N), temp_idx, pv
    integer :: i, j, k, r, c
    
    err = 0
    col_pv = [(i,i=1,N)]
    
    ! Forward elimination
    do i = 1, N-1
      max_val = 0.0d0
      r = i; c = i
      do j = i, N
        do k = i, N
          if (abs(A(j,k)) > max_val) then
            max_val = abs(A(j,k))
            r = j
            c = k
          end if
        end do
      end do
      
      if (max_val < TOL) then
        err = 1
        return
      end if

      ! Swap row
      if (r /= i) then
        temp_vect = A(i,:)
        A(i,:) = A(r,:)
        A(r,:) = temp_vect
        temp = B(i)
        B(i) = B(r)
        B(r) = temp
      end if
      
      ! Swap col
      if (c /= k) then
        temp_vect = A(:,i)
        A(:,i) = A(:,c)
        A(:,c) = temp_vect
        temp_idx = col_pv(i)
        col_pv(i) = col_pv(c)
        col_pv(c) = temp_idx
      end if
      
      do j=k+1,N
        coef = A(j,i)/A(i,i)
        A(j,i:N) = A(j,i:N) - coef * A(i,i:N)
        B(j) = B(j) - coef * B(i)
      end do
    end do

    ! Backward substitution
    if (abs(A(N,N)) < TOL) then
      err = 1
      return
    end if
    
    x(N) = B(N)/A(N,N)
    do i = N-1,1,-1
      x(i) = (B(i) - sum(A(i,i+1:N) * x(i+1:N))) / A(i,i)
    end do

    ! reorder the solution
    do i=1,N
      temp_vect(col_pv(i)) = x(i)
    end do
    x = temp
  end subroutine gauss_eliminate_fpv

end module gauss_elimination
