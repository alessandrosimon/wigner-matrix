module utils
  use ftypes
  use wigner
  use fwigxjpf
  implicit none
    

contains
    function cross(a, b)
        real(dp), dimension(3) :: cross
        real(dp), dimension(3), intent(in) :: a, b

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)
    end function cross

    function get_simulation_box(filename) result(Length)
      character(len=*), intent(in) :: filename
      integer :: fileunit, X, npart
      real(dp) :: Length(3), Lx, Ly, Lz

      open (newunit=fileunit, file=filename, action='read')
      read (unit=fileunit, fmt=*) X, npart, Lx, Ly, Lz
      Length(1) = Lx
      Length(2) = Ly
      Length(3) = Lz
    end function

    ! indexer for g in intermolecular frame
    subroutine get_gim_indexer(lmax, indexer)
      integer, intent(in) :: lmax
      integer, intent(inout), allocatable :: indexer(:,:,:,:,:)
      integer :: idx, l1, l2, chi, n1, n2

      allocate(indexer(0:lmax, 0:lmax, -lmax:lmax, -lmax:lmax, -lmax:lmax))
      idx = 1
      do l1 = 0, lmax
        do l2 = 0, lmax
          do chi = -min(l1, l2), min(l1, l2)
            do n1 = -l1, l1
              do n2 = -l2, l2
                indexer(l1, l2, chi, n1, n2) = idx
                idx = idx+1
              end do
            end do
          end do
        end do
      end do
    end subroutine get_gim_indexer  

    ! indexer for g in lab frame
    subroutine get_g_indexer(lmax, indexer)
      integer, intent(in) :: lmax
      integer, intent(inout), allocatable :: indexer(:,:,:,:,:)
      integer :: idx, l1, l2, l, n1, n2

      allocate(indexer(0:lmax, 0:lmax, -lmax:lmax*2, -lmax:lmax, -lmax:lmax))
      idx = 1
      do l1 = 0, lmax
        do l2 = 0, lmax
          do l = abs(l1-l2), l1+l2
            do n1 = -l1, l1
              do n2 = -l2, l2
                indexer(l1, l2, l, n1, n2) = idx
                idx = idx+1
              end do
            end do
          end do
        end do
      end do
    end subroutine get_g_indexer


    subroutine get_d_indexer(lmax, indexer)
      integer, intent(in) :: lmax
      integer, intent(inout), allocatable :: indexer(:,:,:)
      integer :: l, m, n

      allocate(indexer(0:lmax, -lmax:lmax, -lmax:lmax))
      do l = 0, lmax
        do m = -l, l
          do n = -l, l
            indexer(l, m, n) = Dindexer(l, m, n)
          end do
        end do
      end do

    end subroutine get_d_indexer

    ! transform g_tilde to lab g
    subroutine gt_to_lab(lmax, gt, g)
      integer, intent(in) :: lmax
      complex(dp), intent(inout) :: gt(:, :), g(:, :)
      integer, allocatable, dimension(:,:,:,:,:) :: gindexer, gtindexer
      integer :: l1, l2, l, chi, n1, n2
      real(dp) :: factor

      ! indexer for g and g tilde
      call get_gim_indexer(lmax, gtindexer)
      call get_g_indexer(lmax, gindexer)

      ! r axis is the same
      g(:, 1) = gt(:, 1)

      call fwig_table_init(2*100,9)
      call fwig_temp_init(2*100)


      do l1 = 0, lmax
        do l2 = 0, lmax
          do l = abs(l1-l2), l1+l2
             do chi = -min(l1, l2), min(l1, l2)
                
              factor = fwig3jj(2* l1 , 2* l2 , 2* l , &
                               2* chi, 2* (-chi) , 0)

              ! use luc's convention here
              !              factor = factor * sqrt((2.*l1 + 1.)*(2.*l2 + 1.)) * (2.*l + 1.)
              factor = factor *  (2.*l + 1.)
              
              do n1 = -l1, l1
                do n2 = -l2, l2
                  g(:, 1+gindexer(l1, l2, l, n1, n2)) = g(:,1+gindexer(l1, l2, l, n1, n2)) +&
                  gt(:, 1+gtindexer(l1, l2, chi, n1, n2))*factor
                end do
              end do
            end do
          end do
        end do
      end do

      call fwig_temp_free();
      call fwig_table_free();
      end subroutine gt_to_lab
    


    pure function get_g_dim(lmax) result(dim)
      integer, intent(in) :: lmax
      integer :: dim
      dim = (lmax+1)**2 * (2*lmax+1)**3
    end function get_g_dim

    subroutine init_g_array(garr, nbins, lmax)
      type(ragged_array), allocatable, intent(inout) :: garr(:,:,:,:)
      integer, intent(in) :: lmax, nbins
      integer:: l1, l2, chi, k

      !> maximum number of different indices
      !! more than needed but that's okay
      allocate(garr(nbins, 0:lmax, 0:lmax, -lmax:lmax))

      do l1 = 0, lmax
        do l2 = 0, lmax
          do chi = -min(l1, l2), min(l1, l2)
            do k = 1, nbins
              allocate(garr(k, l1, l2, chi)%m(-l1:l1, -l2:l2))
            end do
          end do
        end do
      end do

    end subroutine

    subroutine linspace(from, to, array)
      real(dp), intent(in) :: from, to
      real(dp), intent(out) :: array(:)
      real(dp) :: range
      integer :: n, i
      n = size(array)
      range = to - from
  
      if (n == 0) return
  
      if (n == 1) then
          array(1) = from
          return
      end if
  
  
      do i=1, n
          array(i) = from + range * (i - 1) / (n - 1)
      end do
  end subroutine
end module utils
