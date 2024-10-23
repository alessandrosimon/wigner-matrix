module wigner
  use ftypes
  implicit none
  integer, allocatable :: ldims(:), cum_dims(:)
  integer :: tot_dim
  !> cache arrays
  type(ragged_array), allocatable::cls(:)
  type(ragged_array), allocatable::dls(:)


contains
  subroutine init_wigner(lmax)
    integer, intent(in) :: lmax
    integer :: i, sum

    !> add one because of D0_00 = unity
    allocate(ldims(lmax+1))
    allocate(cum_dims(lmax+1))


   !> init ldims and sum such that
   !! ldims = 1, 9, 25, ...
   !! sum = \Sigma ldims
    sum = 0
    do i = 0, lmax
       ldims(i+1) = (2*i + 1)**2
       sum = sum + ldims(i+1)
    end do

    !> init cummulative dims:
    !! cum_dims = 1, 10, 35, ...
    cum_dims(1) = 1
    do i = 2, lmax+1
       cum_dims(i) = cum_dims(i-1) + ldims(i)
    end do
    tot_dim = sum

   allocate(cls(lmax))
   allocate(dls(lmax))

   do i = 1, lmax
      ! reuse var 'sum'
      sum = 2*i + 1
      cls(i)%m = get_c(i)
      dls(i)%m = get_d(i)

      ! reverse last dim
      cls(i)%mr = cls(i)%m(:, sum:1:-1)
      dls(i)%mr = dls(i)%m(:, sum:1:-1)


      cls(i)%rmr = cls(i)%m(sum:1:-1, sum:1:-1)
      dls(i)%rmr = dls(i)%m(sum:1:-1, sum:1:-1)


      cls(i)%rm = cls(i)%m(sum:1:-1, :)
      dls(i)%rm = dls(i)%m(sum:1:-1, :)


   end do

  end subroutine init_wigner

  function Dindexer(l, m, mp) result(idx)
   integer, intent(in) :: l, m, mp
   integer :: idx, loffset, moffset, mpoffset

   if (l == 0) then
      idx = 1
      return 
   end if
   if (l == 1) then 
      loffset = 2
   else 
      loffset = cum_dims(l)+1
   end if
   moffset = m+l
   mpoffset = (2*l+1)*(mp+l)
   idx = loffset + moffset + mpoffset   
   
  end function
  
  !> routine splitting the inital rotation matrix
  subroutine get_fg(r, f, g)
    real(dp), intent(in), dimension(0:2, 0:2) :: r
    complex(dp), intent(out), dimension(0:2, 0:2) :: f, g
    real(dp) :: sqrt2

    sqrt2 = sqrt(2.0_dp)

    f(0,:) = [(r(1, 1)+r(0,0))/2.0, r(0,2)/sqrt2, (r(1,1)-r(0,0))/2.]
    f(1,:) = [r(2,0)/sqrt2, R(2,2), -r(2,0)/sqrt2]
    f(2,:) = [(r(1,1)-r(0,0))/2.0, -r(0,2)/sqrt2, (r(1,1)+r(0,0))/2.0]
        
    g(0,:) = [(r(1, 0)-r(0,1))/2.0, r(1,2)/sqrt2, -(r(1,0)+r(0,1))/2.]
    g(1,:) = [-r(2,1)/sqrt2, 0.0_dp, -r(2,1)/sqrt2]
    g(2,:) = [(r(1,0)+r(0,1))/2.0, r(1,2)/sqrt2, (r(0,1)-r(1,0))/2.0]

  end subroutine get_fg

  function get_c(l) result(cl)
    integer, intent(in) :: l
    real(dp), allocatable :: cl(:,:)
    integer :: m,mp

    allocate(cl(-l:l, -l:l))
    do m = -l, l
       do mp = -l, l
          if (abs(m) == l .or. mp == -l .or. mp-1 == -l) then
             cl(m,mp) = 0.0_dp
          else
             cl(m,mp) = sqrt((2.0_dp*(l+m)*(l-m))/((l+mp)*(l+mp-1)))
          endif
       end do
    end do
  end function get_c

  pure function get_dim(l) result(dim)
    integer :: i, dim
    integer, intent(in) :: l
    dim = 1
    do i = 1, l
       dim = dim + (2*i+1)**2
    end do
  end function get_dim

  function get_d(l) result(dl)
    integer, intent(in) :: l
    real(dp), allocatable :: dl(:,:)
    integer :: m,mp

    allocate(dl(-l:l, -l:l))
    do m = -l, l
       do mp = -l, l
          if (m == -l .or. m == (-l+1) .or. mp == -l .or. mp-1==-l) then
             dl(m,mp) = 0.0_dp
          else
             dl(m,mp) = sqrt(1.0_dp*((l+m)*(l+m-1))/((l+mp)*(l+mp-1)))
          endif    
       end do
    end do
  end function get_d

  function wigner_D(lmax, R) result(D_tot)
    integer, intent(in) :: lmax
    integer :: dim, l
    real(dp), intent(in) :: R(3,3)
    complex(dp), allocatable :: D_tot(:)
    !> varables holding the intermediate values of F, G, D
    !! split into the total, left and right part
    complex(dp), allocatable :: D_temp(:,:), F_temp(:,:), G_temp(:,:)
    complex(dp), allocatable :: D_right(:,:), F_right(:,:), G_right(:,:)
    complex(dp), allocatable :: D_left(:,:), F_left(:,:), G_left(:,:)
    !> retain G1 and F1 for later iterations
    complex(dp) :: F1(3,3), G1(3,3)

   ! call init(lmax)
    allocate(D_tot(tot_dim))
    allocate(F_temp(3,3))
    allocate(G_temp(3,3))
    allocate(D_temp(3,3))

    !> split R into F and G
    call get_fg(R, F_temp, G_temp)

    !> construct D1
    D_temp = F_temp + complex(0., 1.0_dp)*G_temp
    F1 = F_temp
    G1 = G_temp
    D_tot = 0.0_dp
    D_tot(1) = 1.0_dp

    if (lmax == 0) then
      return
    end if

    D_tot(2:10) = [D_temp]
    l = 1
    if (lmax == 1) then
       return
    end if


    do while (l < lmax)
       l = l+1
       dim = 2*l+1

       deallocate(D_temp)

       allocate(D_temp(dim,dim))

       allocate(F_right(dim,dim))
       allocate(G_right(dim,dim))
       allocate(D_right(dim,dim))

       allocate(F_left(dim,dim))
       allocate(G_left(dim,dim))
       allocate(D_left(dim,dim))


       F_right = cls(l)%m * padded_h(F_temp, G_temp, [1, 2], [1, 1], [2,0]) &
            + dls(l)%m * padded_h(F_temp, G_temp, [2, 2], [2, 0], [2, 0]) &
            + dls(l)%rm * padded_h(F_temp, G_temp, [0, 2], [0, 2], [2, 0])

       G_right = cls(l)%m * padded_k(F_temp, G_temp, [1, 2], [1, 1], [2, 0]) &
            + dls(l)%m * padded_k(F_temp, G_temp, [2, 2], [2, 0], [2, 0]) &
            + dls(l)%rm * padded_k(F_temp, G_temp, [0, 2], [0, 2], [2, 0])

       F_left = cls(l)%mr * padded_h(F_temp, G_temp, [1, 0], [1, 1], [0, 2]) &
            + dls(l)%mr * padded_h(F_temp, G_temp, [2, 0], [2, 0], [0, 2]) &
            + dls(l)%rmr * padded_h(F_temp, G_temp, [0, 0], [0, 2], [0, 2])

       G_left = cls(l)%mr * padded_k(F_temp, G_temp, [1, 0], [1, 1], [0, 2]) &
            + dls(l)%mr * padded_k(F_temp, G_temp, [2, 0], [2, 0], [0, 2]) &
            + dls(l)%rmr * padded_k(F_temp, G_temp, [0, 0], [0, 2], [0, 2])


       D_right = F_right + complex(0., 1.0_dp) * G_right
       D_left = F_left + complex(0., 1.0_dp) * G_left

       D_temp(:, :l) = D_left(:, :l)
       D_temp(:, l+1:) = D_right(:, l+1:)

       
       deallocate(F_temp, G_temp)
       allocate(F_temp(dim,dim), G_temp(dim, dim))

       F_temp(:, :l) = F_left(:, :l)
       F_temp(:, l+1:) = F_right(:, l+1:)

       G_temp(:, :l) = G_left(:, :l)
       G_temp(:, l+1:) = G_right(:, l+1:)

       D_tot(cum_dims(l)+1:cum_dims(l+1)) = [D_temp]

       deallocate(F_right, G_right, D_right, F_left, G_left, D_left)
    end do
       return
       
  contains
    function padded_h(f, g, idx, pad1, pad2) result(mat)
      !> f and g are both of size 2*(l-1)+1, i.e. those of the earlier iteration
      !! we need to return an array of size 2*l+1 hence we pad f and g in such a way
      !! that the equations from choi99 are fulfilled
      complex(dp), intent(in) :: f(:,:), g(:,:)
      !> the padding variables indicate the number of zeros to be added at front/end of
      !! each dimension (see also numpy.pad)
      !! sum(pad) should always be equal to 2 here (by design)
      integer, intent(in), dimension(2) :: idx, pad1, pad2
      complex(dp), allocatable :: mat(:,:)
      integer :: newdim

      !> f, g are square, first dimension suffices
      newdim = size(f, 1)+2

      !> easier to work with indices starting from zero in this case
      allocate(mat(0:newdim-1, 0:newdim-1))
      mat = 0.0_dp

      mat(pad1(1):newdim-pad1(2) - 1, pad2(1):newdim-pad2(2) - 1) = F1(idx(1)+1, idx(2)+1)*f - G1(idx(1)+1, idx(2)+1)*g
    end function padded_h

    function padded_k(f, g, idx, pad1, pad2) result(mat)
      !> f and g are both of size 2*(l-1)+1, i.e. those of the earlier iteration
      !! we need to return an array of size 2*l+1 hence we pad f and g in such a way
      !! that the equations from choi99 are fulfilled
      complex(dp), intent(in) :: f(:,:), g(:,:)
      !> the padding variables indicate the number of zeros to be added at front/end of
      !! each dimension (see also numpy.pad)
      !! sum(pad) should always be equal to 2 here (by design)
      !! note: we increase idx by one in the code below because F1 and G1
      !! are captured from the outer scope and are indexed starting from 1
      integer, intent(in), dimension(2) :: idx, pad1, pad2
      complex(dp), allocatable :: mat(:,:)
      integer :: newdim

      !> f, g are square, first dimension suffices
      newdim = size(f, 1)+2
      !> easier to work with indices starting from zero in this case
      allocate(mat(0:newdim-1, 0:newdim-1))
      mat = 0.0_dp
      mat(pad1(1):newdim-pad1(2) - 1, pad2(1):newdim-pad2(2) - 1) = F1(idx(1)+1, idx(2)+1)*g + G1(idx(1)+1, idx(2)+1)*f
    end function padded_k
    
  end function wigner_D

  subroutine print_mat(mat)
    implicit none
    complex(dp), intent(in) :: mat(:,:)
    integer :: i
    do i = 1, size(mat, 1)
       write (*, "(*('('sf6.3spf6.2x'i)':x))") mat(i, :)
    end do
  end subroutine print_mat

  subroutine print_matr(mat)
    implicit none
    real(dp), intent(in) :: mat(:,:)
    integer :: i
    do i = 1, size(mat, 1)
       write (*, *) mat(i, :)
    end do
  end subroutine print_matr


end module wigner

! program test
!   use wigner
!   use rotations
!   use utils
!   real(dp), dimension(3, 3) :: a
!   complex(dp), allocatable :: d(:)
!   complex(dp) :: d1(3, 3), d2(5, 5), d3(7, 7), d4(9, 9)
!   integer :: i
  
!   call init_wigner(6)
!   allocate(d(tot_dim))
!   a = matrix_from_euler(0.1_dp, 3.02_dp, 2.23_dp)

!   d = wigner_D(6, a)
!   d1 = reshape(d(2:10), [3, 3])
!   call print_mat(d1)
!   d2 = reshape(d(11:35), [5, 5])
!   !call print_mat(d2)
!   d3 = reshape(d(36:84), [7, 7])
!   call print_mat(d3)
  
!   write (*, *) d(Dindexer(3, 0, 1))
  
  

! end program
