    program add_column
    implicit none
    integer, parameter :: dp = kind(0.d0)  ! Double precision kind
    integer, parameter :: MAX_LINES = 801
    integer, parameter :: MAX_COLUMNS = 6
    real(dp), dimension(MAX_LINES, MAX_COLUMNS) :: data
    character(len=1000) :: line
    integer :: i, j, k, n
    real(dp), dimension(MAX_LINES) :: NG1, NG2
    
    ! Open CSV file
    open(unit=10, file='your_file.csv', status='old', action='read')
    
    ! Skip header line
    read(10, '(a)') line
    
    ! Read data from CSV file
    do i = 1, MAX_LINES
        read(10, *, iostat=n) data(i, 1:MAX_COLUMNS)
        if (n /= 0) exit
    end do
    close(10)
    
    ! Initialize NG1 and NG2
    NG1 = 0
    k = 0
    do i = 1, MAX_LINES
        if (mod(i-1, 4) == 0) then
            k = k + 2
        end if
        NG2(i) = k * 0.1 / 199
    end do
    
    ! Write data to new CSV file
    open(unit=10, file='new_file.csv', status='replace', action='write')
    
    ! Write headers
    write(10, '(7a)') 'Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6', 'NG1', 'NG2'
    
    ! Write data
    do i = 1, MAX_LINES
        write(10, *) (data(i, j), j=1, MAX_COLUMNS), NG1(i), NG2(i)
    end do
    close(10)
    
    end program add_column

