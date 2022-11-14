!subroutine read_parameter
!    use parameter_setting
!    character(len = 20)::filename    ! 文件名称
!    character(len = 80)::msg          ! 打开文件错误时的异常信息
!    integer :: nvals = 0                       ! 读取数据的行数
!    integer :: status1 , status2            ! 打开文件和读取数据的状态信息
!    real(prec) :: value                                   ! 需读取文件中的数据
!
!    filename = 'test.in'
!    write(*,"('The input file name is: 'A)") filename
!    open(UNIT = 3 , FILE = filename , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1 , IOMSG = msg)  ! 打开文件
!    !IOSTAT选项 
!    !   如果OPEN语句执行成功，那么会返回0；
!    !   如果遇到文件结束位置，则返回负整数；
!    !   如果出现错误则返回一个正整数，具体数值视具体的计算机系统而定
!    !IOMSG选项
!    !   OPEN语句执行之后，包含在OPEN语句状态中的字符变量
!    !   文件打开成功，该状态变量的值不会改变
!    !   文件打开失败，变为描述错误的字符信息。
!    
!    openif:if(status1 == 0 ) then                                   ! 文件打开成功
!        readloop:do
!            read(3 , * , IOSTAT = status2) value
!            if(status2 /= 0) exit                                   ! 如果数据读取失败，跳出do循环
!            nvals = nvals + 1
!            write(*,"('Line ',I6,' :Value = ',F20.15)") nvals,value
!        end do readloop
!        
!        readif:if(status2 > 0) then                                             ! 发生了数据读取错误
!                write(*,"('An error occurred reading line ',I6)") nvals+1       ! 交代在哪一行读取数据失败
!            else readif                                                         ! 到达数据结尾
!                write(*,"('end of file reached.There were ' , I6 , ' values in the file.')") nvals  ! 交代已经到达文件结尾
!        end if readif
!        
!    else openif                                                     ! 文件打开失败
!        write(*,"('Error opening file : IOSTAT = ',I6)") status1
!        write(*,*) trim(msg)                                        ! 返回错误信息
!    end if openif
!
!    close(UNIT = 3)     
!
!
!end subroutine read_parameter