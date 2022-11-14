!subroutine read_parameter
!    use parameter_setting
!    character(len = 20)::filename    ! �ļ�����
!    character(len = 80)::msg          ! ���ļ�����ʱ���쳣��Ϣ
!    integer :: nvals = 0                       ! ��ȡ���ݵ�����
!    integer :: status1 , status2            ! ���ļ��Ͷ�ȡ���ݵ�״̬��Ϣ
!    real(prec) :: value                                   ! ���ȡ�ļ��е�����
!
!    filename = 'test.in'
!    write(*,"('The input file name is: 'A)") filename
!    open(UNIT = 3 , FILE = filename , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1 , IOMSG = msg)  ! ���ļ�
!    !IOSTATѡ�� 
!    !   ���OPEN���ִ�гɹ�����ô�᷵��0��
!    !   ��������ļ�����λ�ã��򷵻ظ�������
!    !   ������ִ����򷵻�һ����������������ֵ�Ӿ���ļ����ϵͳ����
!    !IOMSGѡ��
!    !   OPEN���ִ��֮�󣬰�����OPEN���״̬�е��ַ�����
!    !   �ļ��򿪳ɹ�����״̬������ֵ����ı�
!    !   �ļ���ʧ�ܣ���Ϊ����������ַ���Ϣ��
!    
!    openif:if(status1 == 0 ) then                                   ! �ļ��򿪳ɹ�
!        readloop:do
!            read(3 , * , IOSTAT = status2) value
!            if(status2 /= 0) exit                                   ! ������ݶ�ȡʧ�ܣ�����doѭ��
!            nvals = nvals + 1
!            write(*,"('Line ',I6,' :Value = ',F20.15)") nvals,value
!        end do readloop
!        
!        readif:if(status2 > 0) then                                             ! ���������ݶ�ȡ����
!                write(*,"('An error occurred reading line ',I6)") nvals+1       ! ��������һ�ж�ȡ����ʧ��
!            else readif                                                         ! �������ݽ�β
!                write(*,"('end of file reached.There were ' , I6 , ' values in the file.')") nvals  ! �����Ѿ������ļ���β
!        end if readif
!        
!    else openif                                                     ! �ļ���ʧ��
!        write(*,"('Error opening file : IOSTAT = ',I6)") status1
!        write(*,*) trim(msg)                                        ! ���ش�����Ϣ
!    end if openif
!
!    close(UNIT = 3)     
!
!
!end subroutine read_parameter