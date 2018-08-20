module global_parameter
 
  implicit none

  integer, parameter :: BOSONIC_ =  0
  integer, parameter :: FERMIONIC_ = 1

  integer, parameter :: CHANNEL_D = 2
  integer, parameter :: CHANNEL_M = 3
  integer, parameter :: CHANNEL_S = 4
  integer, parameter :: CHANNEL_T = 5

  integer, parameter :: FADDB_    = 6
  integer, parameter :: FADDF_    = 7
  integer, parameter :: FSUBF_    = 8
  integer, parameter :: MINUSF_   = 9
  integer, parameter :: MINUSB_   = 10

  ! mathematic constants
  integer,  parameter :: dp = selected_real_kind(8)
  real(dp), parameter :: pi = acos(-1.d0)
  real(dp), parameter :: Zero = 0.d0, Half = 0.5d0, One = 1.d0, Two = 2.d0  
  complex(dp), parameter :: xi = dcmplx(Zero, One)
  complex(dp), parameter :: Zero_c = dcmplx(Zero, Zero), One_c = dcmplx(One, Zero)
  complex(dp), parameter :: Half_c = dcmplx(Half, Zero), Two_c  = dcmplx(Two, Zero)

  !--- mpi ---
  integer  :: ntasks, master, id, rc
  integer, allocatable :: status(:), send_request(:), recv_request(:)

  ! physical constant
  integer, parameter :: nSpin   = 2
  integer, parameter :: nTau    = 200
  real(dp), parameter :: beta = 2.d0
  real(dp), parameter :: xU   = 4.d0
  real(dp), parameter :: xJ   = 1.d0
  real(dp), parameter :: nParticle = 1.d0

  ! will be determined in Green function
  real(dp) :: mu   = 0.d0

  ! complex(dp), allocatable :: Omega(:)
  ! logical :: Debug_Info = .False.

contains
  !----------------------------------------------------------------------------------------------
  subroutine license_message(num_node)
    implicit none
    integer, intent(in) :: num_node

    write(*, '(a)') ''//achar(27)//'[33m +--------------------------------------------------+'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m         _                                      '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m        (_)      _                              '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m   _   _ _  ____| |_  ___   ____ _   _          '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m  | | | | |/ ___)  _)/ _ \ / ___) | | |         '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m   \ V /| ( (___| |_| |_| | |   | |_| |         '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m    \_/ |_|\____)\___)___/|_|    \__  |         '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m                                (____/          '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[32m    Vienna Computational Tool Depository        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[34m    ==     =             ==           ==        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                                '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m  If you find victory useful in your research,  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m   you should                                   '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m   (i)  cite our papers and the appropriate     '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        references therein AND                  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m   (ii) state in your manuscript/paper that you '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        have used the victory code (or a        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        modified version of it). An appropriate '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        way of acknowledging the use of victory '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        in your publications would be, e.g.     '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m        adding a sentence like                  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                                '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m "The calculation has been performed using the  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m  victory code", followed by the citation to our'//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m  paper: Phys. Rev. B 93, 165103 (2016).        '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                                '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m             Gang Li, Karsten Held @ copyright  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m + '//achar(27)//'[31m                                    2015.06.01  '//achar(27)//'[33m +'//achar(27)//'[0m'
    write(*, '(a)') ''//achar(27)//'[33m +--------------------------------------------------+'//achar(27)//'[0m'
                       
    write(*, *)        
    write(*, '(2x, a, i5, a)') '... MPI starts with', num_node,'  nodes. ...'

  end subroutine license_message
  
  !----------------------------------------------------------------------------------------------  
  subroutine loop_message(ite)
    implicit none
    integer, intent(in) :: ite
    
    ! ---- self-consistent loop monitoring ---
    write(*,*)
    write(*,"(a36)") '!----------------------------------!'
    write(*,"(a24,i5,a7)") '!   Self-Consistent loop', ite, '      !'
    write(*, "(a36)") '!----------------------------------!'
    
  end subroutine loop_message

  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine error_msg(msg)

    character(len=*), intent(in) :: msg

    write(*,*) trim(msg)
    stop
  end subroutine error_msg
  
  !------------------------------------------------------------------------------------------------------ 
  subroutine system_mem_usage
    implicit none

    integer :: getpid
    
    real(kind=8) :: valueRSS

    character(len=200):: filename=' '
    character(len=80) :: line
    character(len=8)  :: pid_char=' '
    integer :: pid, itemp, i
    logical :: ifxst
    
    valueRSS=-1    ! return negative number if not found
    
    !--- get process ID
    
    pid=getpid()
    write(pid_char,'(I8)') pid
    filename='/proc/'//trim(adjustl(pid_char))//'/status'
    
    !--- read system file
    
    inquire (file=filename,exist=ifxst)
    if (.not.ifxst) then
       write (*,*) 'system file does not exist'
       return
    endif
    
    open(unit=100, file=filename, action='read')
    do
       read (100,'(a)',end=120) line
       if (line(1:6).eq.'VmRSS:') then
          read (line(7:),*) itemp
          exit
       endif
    enddo
120 continue

    valueRSS = dble(itemp)
    i = 0
    do while (Int(itemp/1024**i) > 100)
       i = i + 1
    end do
    valueRSS = valueRSS/dble(1024**i)    

    select case (i)
    case (0)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'KB'
    case (1)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'MB'
    case (2)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'GB'
    case (3)
      write(*, '(a,f8.4, a3)'), 'Memory used on each node is roughly:', valueRSS, 'TB'
    case default
      write(*, '(a)') 'Memory used more than 1024 TB'
    end select
    close(100)

  end subroutine system_mem_usage

end module global_parameter
