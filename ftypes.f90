module ftypes
    use iso_fortran_env, only: sp=>real32, dp=>real64, qp=>real128
    use iso_fortran_env, only: int8, int16, int32, int64

    !> define array of matrices with different size to later cache the values of 
    !! the matrices c(l) and d(l)
    type ragged_array  
        complex(dp), allocatable::m(:, :)
        complex(dp), allocatable::rm(:, :)
        complex(dp), allocatable::mr(:, :)
        complex(dp), allocatable::rmr(:, :)
    end type ragged_array


end module ftypes   