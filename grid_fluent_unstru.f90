subroutine fluent_unstru_input
     
    !-----------------------------------------------------------------------------
    !
    !   ����������fluent *.cas��ʽ�ļ�
    !   ������*.cas��ʽ�ļ�˵����Fluent_File_Doc.txt
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none 

    call read_mesh          !----�����ȡ
    call get_stru           !----���ݽṹ����
    call solve_SPs_local    !----����Jacobi�����ȫ������ֵ  
    call solve_FPs_local    !----ͨ����Jacobi�����ȫ������ֵ

end subroutine fluent_unstru_input
    
subroutine read_mesh
    use real_precision
    use type_module
    use global_var
    use parameter_setting,only:mesh_file
    implicit none 
    
    integer :: IOfile = 10000               !----�����ļ�
    real(prec) :: coor_x,coor_y             !----��������
    
    character(len = 100):: filename 
    character(len = 1 ) :: char_temp                    
    character(len = 100):: char_temp1,char_temp_1ast
    character(len = 100) :: char_line_1,char_line_2,char_line_sub(20)
    character(len = 16 ) :: bc
    character(len = 10 ):: char_temp2
    character(len = 8)  :: char_index(4)    !----������飬16����
    integer :: int_index(4)                 !----������飬10����
    
    integer nnodes1,nsides1,ncells1,bdside1
    integer len_char
    integer i,j,k,l,istart,iend,nsize
    
    logical alive                           !----�ļ�������Ϣ
    character(len = 100)::msg               !----���ļ�����ʱ���쳣��Ϣ
    integer :: status1 , status2=0          !----���ļ��Ͷ�ȡ���ݵ�״̬��Ϣ

    integer node1,node2,cell1,cell2,temp
    real(prec) :: coor_vertex(4,2)
     
    
    !-----------------------------------------------------------------------------
    !
    !   �������������ļ�
    !
    !-----------------------------------------------------------------------------
    
    write(*,*) 'Mesh Filename:',mesh_file
    filename = '.\mesh\'//mesh_file
    
    inquire(file=filename,exist=alive)
    if(alive)then
        open(UNIT = IOfile , FILE = filename , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1 , IOMSG = msg)
        if(status1 == 0 ) then              !----�ļ��򿪳ɹ�
            do i = 1,6
                read(IOfile,*)              !----����ǰ6��
            end do
        else                                !----�ļ���ʧ��
            write(*,"('Error opening file : IOSTAT = ',I6)") status1
            write(*,*) trim(msg)            !----���ش�����Ϣ
            write(*,*)'0000'
            stop
        end if     
    else
        write(*,*)'Error! No mesh file exist, occur in fluent_unstru_input. '
        stop
    end if
      
    !-----------------------------------------------------------------------------
    !
    !   ��������ȡnode side cell����Ŀ
    !   ������*.cas�ļ���ͷ���ݴ�����ͬ�������ù̶���ʽ�ж�.�������ݣ����ݱ���жϣ��ο�Fluent_File_Doc.txt
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
     
    ! ��ȡ�ܵ���nnodes
    read(IOfile,*) char_temp,char_temp1
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(19:len_char)
      
    read(char_temp_1ast,*) nnodes1
 
    ! ��ȡ�ܱ���nsides
    do i = 1, 2
        read(IOfile,*)
    end do
    read(IOfile,*) char_temp,char_temp1   
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(25:len_char)              
    read(char_temp_1ast,*) nsides1

    ! ��ȡ�ܱ߽����bdside
    read(IOfile,*) char_temp,char_temp1
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(25:len_char) 
    read(char_temp_1ast,*) bdside1
    
    ! ��ȡ�ܵ�Ԫ��nside
    do i = 1,3
       read(IOfile,*)
    end do
    read(IOfile,*) char_temp,char_temp1
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(24:len_char)           
    read(char_temp_1ast,*) ncells1
    
    ! ���ݸ�ȫ�ֱ���
    nnodes = nnodes1
    nsides = nsides1
    ncells = ncells1
    nbdsides = bdside1
    write(*,"('   nodes:'I8,'  sides:'I8,'  cells:'I8,'  nbdsides:'I8)")nnodes,nsides,ncells,nbdsides
       
    !-----------------------------------------------------------------------------
    !
    !   ���������꼰��Ŷ�ȡ
    !
    !-----------------------------------------------------------------------------
    
    ! ����ṹ�������ڴ沢��ʼ��
    call allo_stru_record   
    allocate(xy_coor(nnodes+nbdsides*2,2))
    
    xy_coor = 0.0_prec
    do i = 1,ncells
        cellset(i).nodes = 0
        cellset(i).sides = 0
        cellset(i).nearcells = 0
    end do
    
    ! ���������coord
    do i = 1,6
       read(IOfile,*)
    end do
    do i = 1,nnodes1
        read(IOfile,*)coor_x,coor_y
        xy_coor(i,1) = coor_x
        xy_coor(i,2) = coor_y
        
        !�Զ�����������Χ��ʡȥ�ֶ���������������ʺ��ı��μ�������
        if(i==1)then
            xl = coor_x
            xr = coor_x
            yl = coor_y
            yr = coor_y      
        end if
        xl = min(xl,coor_x)
        xr = max(xr,coor_x)
        yl = min(yl,coor_y)
        yr = max(yr,coor_y)
    end do
    xlong = xr - xl
    ylong = yr - yl
    write(*,*) "Calculation area ",xl,xr,yl,yr
         
    !-----------------------------------------------------------------------------
    !
    !   ��������ÿ���ߵı������
    !
    !-----------------------------------------------------------------------------
     
    char_index = '0'

    do while(status2 == 0 )
        READ(IOfile , * , IOSTAT = status2) char_index(1:2)                 
        if(trim(char_index(1)) == '(13'.and.char_index(3)=='0')then         !----13��ʾ�ߣ���Fluent_File_Doc.txt
            
            ! �˺����У����߽�����
            Backspace(IOfile)
            Backspace(IOfile)
            read(IOfile ,* , IOSTAT = status2)char_line_1, char_line_2
            call StringSplit(char_line_2,' ',char_line_sub,nsize)           !----�ַ����ָ�
            bc = char_line_sub(nsize)

            ! ��ȡ�ߵ���Ϣ����������㣬�յ�
            READ(IOfile , * , IOSTAT = status2) char_index    
            read(char_index(3),"(Z8)")istart                                !----16����ת��Ϊ10����
            read(char_index(4),"(Z8)")iend 

            do i = istart,iend
                sideset(i).bc = bc                                          !----ÿ���ߵı߽��������ڲ�����Ϊ '0'
                READ(IOfile , * , IOSTAT = status2) char_index              !----����ÿ���ߵı������
                READ(char_index , "(Z8)" , IOSTAT = status2) int_index(:)   !----(��㣬�յ㣬�ҵ�Ԫ����Ԫ)
                
                ! �ѵ�ͱ߷��䵽��Ԫ�ṹ��ͱ߽ṹ����
                node1 = int_index(1)
                node2 = int_index(2)
                cell1 = int_index(3)
                cell2 = int_index(4)
                
                sideset(i).nodes = int_index(1:2)
                sideset(i).nearcells = int_index(3:4)
                
                ! ��ⵥԪ�ṹ�����Ƿ��Ѵ��ڵ���
                if( cell1/=0)then
                    check: do j = 1,4                   
                        if(cellset(cell1).nodes(j) == 0)then
                            if(j == 1)then
                                !�ȴ���1,2��
                                cellset(cell1).nodes(1) = node1
                                cellset(cell1).nodes(2) = node2
                                exit check
                            else
                                if(node1 == cellset(cell1).nodes(1))then
                                    cellset(cell1).nodes(4) = node2
                                elseif(node1 == cellset(cell1).nodes(2))then
                                    cellset(cell1).nodes(3) = node2
                                elseif(node2 == cellset(cell1).nodes(1))then
                                    cellset(cell1).nodes(4) = node1
                                elseif(node2 == cellset(cell1).nodes(2))then
                                    cellset(cell1).nodes(3) = node1
                                else
                                    cellset(cell1).nodes(3) = node1
                                    cellset(cell1).nodes(4) = node2
                                end if
                                exit check
                            end if
                        end if                    
                    end do check
                end if
                
                if(cell2/=0)then
                    check2: do j = 1,4               
                        if(cellset(cell2).nodes(j) == 0  )then
                            if(j == 1)then
                                cellset(cell2).nodes(1) = node1
                                cellset(cell2).nodes(2) = node2
                                exit check2
                            else
                                if(node1 == cellset(cell2).nodes(1))then
                                    cellset(cell2).nodes(4) = node2
                                elseif(node1 == cellset(cell2).nodes(2))then
                                    cellset(cell2).nodes(3) = node2
                                elseif(node2 == cellset(cell2).nodes(1))then
                                    cellset(cell2).nodes(4) = node1
                                elseif(node2 == cellset(cell2).nodes(2))then
                                    cellset(cell2).nodes(3) = node1                 
                                else
                                    cellset(cell2).nodes(3) = node1
                                    cellset(cell2).nodes(4) = node2
                                end if
                                exit check2
                            end if
                        endif         
                    end do check2                      
                end if   
            end do      
            
            char_index = '0'   !----����Backspace(IOfile)������������ѭ��     

        end if
    end do

    close(IOfile)
   
end subroutine read_mesh

subroutine get_stru
 
    !-----------------------------------------------------------------------------
    !
    !   ���������ṹ������
    !   �����������Ե�Ԫ�����������š���ߡ��ڵ�Ԫ
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use type_module
    use global_var
    implicit none 
   
    integer i,j,k,l,istart,iend
    integer indexCellL,indexCellR
    integer node1,node2,cell1,cell2,temp
    real(prec) :: coor_vertex(4,2)
         
    !-----------------------------------------------------------------------------
    !
    !   ��������cellset��nodes����������ʱ��
    !   ������������һ�����صı�ţ���ȡ������ʱ�ȽϷ��㡣һ������Ϊ1����ʱ����ת1-2-3-4-1  
    !   ������  �ж����Ҳ࣬�������
    !           P��Q<0˵��P��Q����ʱ�뷽λ
    !           P��Q>0˵��P��Q��˳ʱ�뷽λ
    !           P��Q=0˵��P��Q����(ͬ�򡢷���)
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
     
    do i = 1,ncells
        do j =1,4
            coor_vertex(j,:) = xy_coor(cellset(i).nodes(j),:)
        end do
        
        ! ��ֵ����
        call sort_qua(coor_vertex(:,:))
      
        ! ��ֵ�������,��ʱ��
        do j = 1,4
            do k = 1,4
                if(coor_vertex(j,1)==xy_coor(cellset(i).nodes(k),1).and.coor_vertex(j,2)==xy_coor(cellset(i).nodes(k),2))then
                    temp = cellset(i).nodes(j)   
                    cellset(i).nodes(j) = cellset(i).nodes(k) 
                    cellset(i).nodes(k) = temp
                end if
            end do
        end do      
    end do
    
    ! �����ߺ��ڵ�Ԫ
    cycle1:do i = 1,nsides   
        indexCellL = sideset(i).nearcells(1)
        indexCellR = sideset(i).nearcells(2)
        if(indexCellL /= 0)then
            cycle2:do j = 1,4
                if(sideset(i).nodes(1)==cellset(indexCellL).nodes(j))then
                    cycle3:do k = 1,4                   
                        if(sideset(i).nodes(2)==cellset(indexCellL).nodes(k))then
                            if((j==1.and.k==4).or.(j==4.and.k==1))then
                                cellset(indexCellL).sides(4)=i
                                cellset(indexCellL).nearcells(4)=indexCellR
                            else
                                cellset(indexCellL).sides(min(j,k))=i
                                cellset(indexCellL).nearcells(min(j,k))=indexCellR
                            end if
                            exit cycle2
                        end if
                    end do cycle3
                end if  
            end do cycle2
   
        end if
        
        if(indexCellR/=0)then
            cycle22:do j = 1,4
                if(sideset(i).nodes(1)==cellset(indexCellR).nodes(j))then
                    cycle33:do k = 1,4                   
                        if(sideset(i).nodes(2)==cellset(indexCellR).nodes(k))then
                            if((j==1.and.k==4).or.(j==4.and.k==1))then
                                cellset(indexCellR).sides(4)=i
                                cellset(indexCellR).nearcells(4)=indexCellL
                            else
                                cellset(indexCellR).sides(min(j,k))=i
                                cellset(indexCellR).nearcells(min(j,k))=indexCellL
                            end if
                            exit cycle22
                        end if
                    end do cycle33
                end if 
            end do cycle22

        end if

    end do cycle1
    
    ! index
    do i = 1,ncells
        cellset(i).index = i
    end do
    
end subroutine get_stru

subroutine sort_qua(arr)
     
    !-----------------------------------------------------------------------------
    !
    !   �������Ե�Ԫ���������ֵ��������
    !   �������Ե�Ԫ�ĸ�����б�ţ�һ�ֲ��ԣ���Ԫ���½Ƕ���Ϊ1��Ȼ����ʱ����
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none
    
    real(prec) :: arr(4,2),arr_temp(4,2),tmp(2)
    integer i,j,k
    real(prec),dimension(4,2) :: ni
    real(prec),external :: Cross_product
    
    ! ��xֵ���� С->��
    call Bubble_Sort(arr(:,:),4)

    ! ��xֵ��С��������ѡ��yֵ��С����Ϊһ����Ԫ����ʼ��
    if(arr(1,2)>arr(2,2))then
        tmp=arr(1,:)
        arr(1,:)=arr(2,:)
        arr(2,:)=tmp
    end if
    
    ! ��yֵ��С��������ѡ��xֵ��С����Ϊһ����Ԫ����ʼ��
    !if(arr(1,1)>arr(2,1))then
    !    tmp=arr(1,:)
    !    arr(1,:)=arr(2,:)
    !    arr(2,:)=tmp
    !end if
    

    !-----------------------------------------------------------------------------
    !
    !   ������������
    !   ����������12 �� ����14 ��ˡ�����������ϵ�ҵ��Խǵ㡣����1ȷ�������μ����ǶԽǵ�
    !   ���ߣ�gqShi 
    !   ��ʷ��** //2021.12.08
    !
    !-----------------------------------------------------------------------------
     
    do i = 1,4
        ni(i,:) = arr(i,:)-arr(1,:)
    end do
    arr_temp(1,:) = arr(1,:)
    cycle1:do i = 2,4
        cycle2:do j = 2,4
            if(j/=i)then
                cycle3:do k =2,4
                    if(k/=j)then
                        if(Cross_product(ni(i,:) ,ni(j,:) )>0 .and. Cross_product(ni(i,:) ,ni(k,:))>0)then
                            if(Cross_product(ni(j,:),ni(k,:))>0)then
                                arr_temp(2,:) = arr(i,:)
                                arr_temp(3,:) = arr(j,:)
                                arr_temp(4,:) = arr(k,:)   
                            else
                                arr_temp(2,:) = arr(i,:)
                                arr_temp(3,:) = arr(k,:)
                                arr_temp(4,:) = arr(j,:)                    
                            end if
                            exit cycle1                     
                        end if
                    end if
                end do cycle3            
            end if
        end do cycle2
    end do cycle1
    arr = arr_temp
    
end subroutine sort_qua

function Cross_product(n1,n2)

    !���
    use real_precision
    implicit none
    
    real(prec) :: Cross_product
    real(prec) ::n1(2),n2(2)
    
    Cross_product = n1(1)*n2(2) - n1(2)*n2(1)
    return
end function Cross_product

subroutine Bubble_Sort(arr,N)

    !ð������
    use real_precision
    implicit none
    
    integer :: i,j,N
    real(prec) :: tmp(2)
    real(prec) :: arr(N,2)
      
    do i=1,N-1
        do j=i+1,N
            !!��xֵ
            if(arr(j,1)<arr(i,1))then
                tmp=arr(i,:)
                arr(i,:)=arr(j,:)
                arr(j,:)=tmp
            end if
            !!��yֵ
            !if(arr(j,2)<arr(i,2))then
            !    tmp=arr(i,:)
            !    arr(i,:)=arr(j,:)
            !    arr(j,:)=tmp
            !end if
        end do
    end do
end subroutine Bubble_Sort