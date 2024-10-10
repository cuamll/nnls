program main
  use iso_c_binding
  use nnls
  implicit none
  integer(kind=CI) :: i, j, k, n, m, maxdim, mode, maxiter, failures, outunit
  integer(kind=c_long) ::  start_time, end_time, count_rate
  real(kind=CF) :: r, res, tol, diff, xdiff
  real(kind=CF), dimension(:), allocatable :: b, x, x_ref
  real(kind=CF), dimension(:, :), allocatable :: a
  character(len=100) :: outfile

  failures = 0_CI
  maxdim = 30_CI
  call random_init(.true., .true.)
  call system_clock(start_time, count_rate)
  do i = 1, 1000
    call random_number(r)
    ! ensure matrix size is between 2 x 2 and maxdim x maxdim
    n = 2_CI + floor((maxdim - 2_CI) * r) 
    call random_number(r)
    m = 2_CI + floor((maxdim - 2_CI) * r) 
    mode = 0_CI
    res = 0.0_CF
    maxiter = 3 * n
    tol = 1.0e-6_CF

    allocate(b(m), source=0.0_CF)
    allocate(x(n), source=0.0_CF)
    allocate(x_ref(n))
    allocate(a(m, n))
    call random_number(a)
    call random_number(x_ref)
    b = matmul(a, x_ref)
    call solve(a, b, x, m, n, mode, res, maxiter, tol)
    diff = norm2(matmul(A, x) - b)
    xdiff = norm2(x_ref - x)
    if ((mode.eq.-1_CI).or.(diff.gt.tol)) then
      write(outfile, '(a, i0.4, a)') "out/info_", i, ".txt"
      open(newunit=outunit, file=outfile)
      failures = failures + 1
      write(outunit, *) "iteration ", i
      write(outunit, *) "x_ref = ", x_ref
      write(outunit, *) "x = ", x
      write(outunit, *) "shape(A) = ", shape(A)
      write(outunit, *) "A = "
      do j = 1, m
        write(outunit, *) (A(j, k), k = 1, n)
      end do
      write(outunit, *)
      write(outunit, *) "Ax = ", matmul(A, x)
      write(outunit, *) "b = ", b
      write(outunit, *) "diff = ", diff
      close(outunit)
    else
      write(*, '(a, i4, a, G10.3, a, G10.3)') "i = ", i,&
        ": |Ax - b| = ", diff, ", |x_ref - x| = ", xdiff
    end if
    deallocate(b)
    deallocate(x)
    deallocate(x_ref)
    deallocate(a)
  end do
  call system_clock(end_time)
  write(*, *) "number of failures: ", failures
  write(*, '(a, f10.6)') "time elapsed (s) = ",&
    real(end_time - start_time) / real(count_rate)

end program main
