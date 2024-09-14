      subroutine ran(r,inseed)
      integer values(1:8), k
      integer, dimension(:), allocatable :: seed
      real*8 r
      integer inseed

      call date_and_time(values=values)

      call random_seed(size=k)
      allocate(seed(1:k))
      if (inseed.eq.0) then
      seed(:) = values(8)
      else
      seed(:) = inseed
      endif
      call random_seed(put=seed)
      call random_number(r)

      return
      END
