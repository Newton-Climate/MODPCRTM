      program cat
      implicit none
      integer n
      character*10 string, number
      string = 'Berkeley'
      n = 102004
      string = trim(string)
      write(number, '(i0)') n
      write(*,*), trim(string)//trim(number)
      end program cat
