! 13/02/2024
! @author: callum
! refs: 
! 1. scipy.optimize.nnls
! 2. Lawson & Hanson 
!    doi.org/10.1137/1.9781611971217
! 3. Bro & de Jong
!    doi.org/10.1002/(SICI)1099-128X(199709/10)11:5<393::AID-CEM483>3.0.CO;2-L
module nnls
  use iso_fortran_env
  use iso_c_binding
  implicit none
  public :: update_s, solve
  integer, parameter :: CF = c_double
  integer, parameter :: CI = c_int
  integer, parameter :: CB = c_bool

  contains

    subroutine update_s(s, ata, atb, pr, n, n_true, s_p_min)&
        bind(C, name='update_s')
      implicit none
      integer(kind=CI) :: n, n_true, i, info, lwork
      real(kind=CF), dimension(n, n) :: ata
      real(kind=CF), dimension(n) :: atb, s
      logical(kind=CB), dimension(n) :: pr
      real(kind=CF) :: s_p_min
      integer(kind=CI), dimension(n_true) :: p_i, ipiv
      real(kind=CF), dimension(n_true, n_true) :: atap
      real(kind=CF), dimension(n_true) :: atbp
      real(kind=CF), dimension(:), allocatable :: work

      if (size(ata, 1).ne.size(ata, 2)) then
        write (*,*) "ata not square"
        stop
      end if
      if (size(ata, 1).ne.size(pr)) then
        write (*,*) "ata and p not same size"
        stop
      end if

      ! indices where p is true
      p_i = pack([(i, i = 1_CI, size(pr))], pr)
      atap = ata(p_i, p_i)
      atbp = atb(p_i)

      allocate(work(1))
      call dsysv("U", n_true, 1, atap, n_true,&
                 ipiv, atbp, n_true, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      call dsysv("U", n_true, 1, atap, n_true,&
                 ipiv, atbp, n_true, work, lwork, info)
      ! atbp is now s[p]
      s_p_min = huge(0.0_CF)
      do i = 1, n_true
        s(p_i(i)) = atbp(i)
        if (s(p_i(i)).lt.s_p_min) then
          s_p_min = s(p_i(i))
        end if
      end do
      deallocate(work)

    end subroutine update_s

    subroutine solve(A, b, x, m, n, mode, res, maxiter, tol)&
        bind(C, name='solve')
      implicit none
      integer(kind=CI) :: m, n, iter, i, mode, k, maxiter
      real(kind=CF) :: tol, s_p_min, alpha, alpha_min, res
      real(kind=CF), dimension(m, n) :: A
      real(kind=CF), dimension(m) :: b
      real(kind=CF), dimension(n) :: x
      real(kind=CF), dimension(:, :), allocatable :: ata
      real(kind=CF), dimension(:), allocatable :: atb, resid, s
      logical(kind=CB), dimension(:), allocatable :: pr

      if (size(b).ne.m) then
        write(*, *) size(A, 1), size(A, 2), size(b)
        write(*, *) "A and b dimensions incompatible"
        stop
      end if
      if (size(x).ne.n) then
        write(*, *) size(A, 1), size(A, 2), size(x)
        write(*, *) "A and x dimensions incompatible"
        stop
      end if

      allocate(ata(n, n), source=0.0_CF)
      ata = matmul(transpose(A), A)
      allocate(atb(n), source=0.0_CF)
      allocate(resid(n), source=0.0_CF)
      atb = matmul(b, A)
      resid = atb
      x = 0.0_CF
      allocate(pr(n))
      pr = .false._CB
      allocate(s(n), source=0.0_CF)

      mode = 1
      iter = 0
      do while ((.not.all(pr)).and.&
        (any(merge(resid, 0.0_CF,&
        (pr.eqv..false._CB)).gt.tol)))

        where (pr) resid = -huge(0.0_CF)
        k = maxloc(resid, 1, kind=CI)
        pr(k) = .true._CB ! must specify dim to get scalar

        s = 0.0_CF
        call update_s(s, ata, atb, pr, n, count(pr), s_p_min)

        do while ((iter.lt.maxiter).and.(s_p_min.le.tol))

          alpha_min = huge(0.0_CF)
          do i = 1, n
            if ((pr(i)).and.s(i).le.tol) then
              alpha = (x(i) / (x(i) - s(i))) 
              if (alpha.lt.alpha_min) then
                alpha_min = alpha
              end if
            end if
          end do
          x = x * (1.0_CF - alpha_min)
          x = x + (alpha_min * s)
          where (x.lt.tol) pr = .false._CB
          call update_s(s, ata, atb, pr, n, count(pr), s_p_min)
          where (.not.pr) s = 0.0_CF
          iter = iter + 1_CI
          
        end do

        x = s
        resid = atb - matmul(ata, x)
        if (iter.ge.maxiter) then
          x = 0.0_CF
          res = 0.0_CF
          mode = -1_CI
          deallocate(ata)
          deallocate(atb)
          deallocate(resid)
          deallocate(pr)
          deallocate(s)
          return
        end if

      end do

      res = norm2(matmul(A, x) - b)
      deallocate(ata)
      deallocate(atb)
      deallocate(resid)
      deallocate(pr)
      deallocate(s)

    end subroutine solve

end module nnls
