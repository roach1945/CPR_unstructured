subroutine fluent_unstru_input
     
    !-----------------------------------------------------------------------------
    !
    !   方法：导入fluent *.cas格式文件
    !   描述：*.cas格式文件说明见Fluent_File_Doc.txt
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none 

    call read_mesh          !----网格读取
    call get_stru           !----数据结构调整
    call solve_SPs_local    !----求解点Jacobi矩阵和全局坐标值  
    call solve_FPs_local    !----通量点Jacobi矩阵和全局坐标值

end subroutine fluent_unstru_input
    
subroutine read_mesh
    use real_precision
    use type_module
    use global_var
    use parameter_setting,only:mesh_file
    implicit none 
    
    integer :: IOfile = 10000               !----输入文件
    real(prec) :: coor_x,coor_y             !----横纵坐标
    
    character(len = 100):: filename 
    character(len = 1 ) :: char_temp                    
    character(len = 100):: char_temp1,char_temp_1ast
    character(len = 100) :: char_line_1,char_line_2,char_line_sub(20)
    character(len = 16 ) :: bc
    character(len = 10 ):: char_temp2
    character(len = 8)  :: char_index(4)    !----编号数组，16进制
    integer :: int_index(4)                 !----编号数组，10进制
    
    integer nnodes1,nsides1,ncells1,bdside1
    integer len_char
    integer i,j,k,l,istart,iend,nsize
    
    logical alive                           !----文件存在信息
    character(len = 100)::msg               !----打开文件错误时的异常信息
    integer :: status1 , status2=0          !----打开文件和读取数据的状态信息

    integer node1,node2,cell1,cell2,temp
    real(prec) :: coor_vertex(4,2)
     
    
    !-----------------------------------------------------------------------------
    !
    !   方法：打开网格文件
    !
    !-----------------------------------------------------------------------------
    
    write(*,*) 'Mesh Filename:',mesh_file
    filename = '.\mesh\'//mesh_file
    
    inquire(file=filename,exist=alive)
    if(alive)then
        open(UNIT = IOfile , FILE = filename , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1 , IOMSG = msg)
        if(status1 == 0 ) then              !----文件打开成功
            do i = 1,6
                read(IOfile,*)              !----跳过前6行
            end do
        else                                !----文件打开失败
            write(*,"('Error opening file : IOSTAT = ',I6)") status1
            write(*,*) trim(msg)            !----返回错误信息
            write(*,*)'0000'
            stop
        end if     
    else
        write(*,*)'Error! No mesh file exist, occur in fluent_unstru_input. '
        stop
    end if
      
    !-----------------------------------------------------------------------------
    !
    !   方法：读取node side cell的数目
    !   描述：*.cas文件开头内容大致相同，可以用固定句式判断.后面内容，根据标号判断，参考Fluent_File_Doc.txt
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
     
    ! 读取总点数nnodes
    read(IOfile,*) char_temp,char_temp1
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(19:len_char)
      
    read(char_temp_1ast,*) nnodes1
 
    ! 读取总边数nsides
    do i = 1, 2
        read(IOfile,*)
    end do
    read(IOfile,*) char_temp,char_temp1   
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(25:len_char)              
    read(char_temp_1ast,*) nsides1

    ! 读取总边界边数bdside
    read(IOfile,*) char_temp,char_temp1
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(25:len_char) 
    read(char_temp_1ast,*) bdside1
    
    ! 读取总单元数nside
    do i = 1,3
       read(IOfile,*)
    end do
    read(IOfile,*) char_temp,char_temp1
    len_char = len_trim(char_temp1)
    char_temp_1ast = char_temp1(24:len_char)           
    read(char_temp_1ast,*) ncells1
    
    ! 传递给全局变量
    nnodes = nnodes1
    nsides = nsides1
    ncells = ncells1
    nbdsides = bdside1
    write(*,"('   nodes:'I8,'  sides:'I8,'  cells:'I8,'  nbdsides:'I8)")nnodes,nsides,ncells,nbdsides
       
    !-----------------------------------------------------------------------------
    !
    !   方法：坐标及编号读取
    !
    !-----------------------------------------------------------------------------
    
    ! 分配结构体数组内存并初始化
    call allo_stru_record   
    allocate(xy_coor(nnodes+nbdsides*2,2))
    
    xy_coor = 0.0_prec
    do i = 1,ncells
        cellset(i).nodes = 0
        cellset(i).sides = 0
        cellset(i).nearcells = 0
    end do
    
    ! 读点的坐标coord
    do i = 1,6
       read(IOfile,*)
    end do
    do i = 1,nnodes1
        read(IOfile,*)coor_x,coor_y
        xy_coor(i,1) = coor_x
        xy_coor(i,2) = coor_y
        
        !自动检索计算域范围，省去手动更改输入参数。适合四边形计算区域
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
    !   方法：读每条边的编号内容
    !
    !-----------------------------------------------------------------------------
     
    char_index = '0'

    do while(status2 == 0 )
        READ(IOfile , * , IOSTAT = status2) char_index(1:2)                 
        if(trim(char_index(1)) == '(13'.and.char_index(3)=='0')then         !----13表示边，见Fluent_File_Doc.txt
            
            ! 退后两行，读边界条件
            Backspace(IOfile)
            Backspace(IOfile)
            read(IOfile ,* , IOSTAT = status2)char_line_1, char_line_2
            call StringSplit(char_line_2,' ',char_line_sub,nsize)           !----字符串分割
            bc = char_line_sub(nsize)

            ! 读取边的信息。边数，起点，终点
            READ(IOfile , * , IOSTAT = status2) char_index    
            read(char_index(3),"(Z8)")istart                                !----16进制转化为10进制
            read(char_index(4),"(Z8)")iend 

            do i = istart,iend
                sideset(i).bc = bc                                          !----每条边的边界条件，内部边设为 '0'
                READ(IOfile , * , IOSTAT = status2) char_index              !----关于每条边的编号内容
                READ(char_index , "(Z8)" , IOSTAT = status2) int_index(:)   !----(起点，终点，右单元，左单元)
                
                ! 把点和边分配到单元结构体和边结构体中
                node1 = int_index(1)
                node2 = int_index(2)
                cell1 = int_index(3)
                cell2 = int_index(4)
                
                sideset(i).nodes = int_index(1:2)
                sideset(i).nearcells = int_index(3:4)
                
                ! 检测单元结构体中是否已存在点或边
                if( cell1/=0)then
                    check: do j = 1,4                   
                        if(cellset(cell1).nodes(j) == 0)then
                            if(j == 1)then
                                !先存入1,2点
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
            
            char_index = '0'   !----控制Backspace(IOfile)，以免陷入死循环     

        end if
    end do

    close(IOfile)
   
end subroutine read_mesh

subroutine get_stru
 
    !-----------------------------------------------------------------------------
    !
    !   方法：求解结构体属性
    !   描述：包括对单元索引、顶点编号、侧边、邻单元
    !   作者：gqShi 
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
    !   方法：对cellset的nodes重新排序，逆时针
    !   描述：给顶点一个当地的编号，在取出数据时比较方便。一般左下为1，逆时针旋转1-2-3-4-1  
    !   方法：  判断左右侧，向量叉乘
    !           P×Q<0说明P在Q的逆时针方位
    !           P×Q>0说明P在Q的顺时针方位
    !           P×Q=0说明P、Q共线(同向、反向)
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
     
    do i = 1,ncells
        do j =1,4
            coor_vertex(j,:) = xy_coor(cellset(i).nodes(j),:)
        end do
        
        ! 对值排序
        call sort_qua(coor_vertex(:,:))
      
        ! 按值调整编号,逆时针
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
    
    ! 分配侧边和邻单元
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
    !   方法：对单元顶点的坐标值进行排序
    !   描述：对单元四个点进行编号，一种策略，单元左下角顶点为1，然后逆时针编号
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none
    
    real(prec) :: arr(4,2),arr_temp(4,2),tmp(2)
    integer i,j,k
    real(prec),dimension(4,2) :: ni
    real(prec),external :: Cross_product
    
    ! 按x值排序 小->大
    call Bubble_Sort(arr(:,:),4)

    ! 从x值较小的两组数选择y值最小的作为一个单元的起始点
    if(arr(1,2)>arr(2,2))then
        tmp=arr(1,:)
        arr(1,:)=arr(2,:)
        arr(2,:)=tmp
    end if
    
    ! 从y值较小的两组数选择x值最小的作为一个单元的起始点
    !if(arr(1,1)>arr(2,1))then
    !    tmp=arr(1,:)
    !    arr(1,:)=arr(2,:)
    !    arr(2,:)=tmp
    !end if
    

    !-----------------------------------------------------------------------------
    !
    !   方法：顶点编号
    !   描述：向量12 与 向量14 叉乘。利用向量关系找到对角点。顶点1确定，依次假设是对角点
    !   作者：gqShi 
    !   历史：** //2021.12.08
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

    !叉乘
    use real_precision
    implicit none
    
    real(prec) :: Cross_product
    real(prec) ::n1(2),n2(2)
    
    Cross_product = n1(1)*n2(2) - n1(2)*n2(1)
    return
end function Cross_product

subroutine Bubble_Sort(arr,N)

    !冒泡排序
    use real_precision
    implicit none
    
    integer :: i,j,N
    real(prec) :: tmp(2)
    real(prec) :: arr(N,2)
      
    do i=1,N-1
        do j=i+1,N
            !!按x值
            if(arr(j,1)<arr(i,1))then
                tmp=arr(i,:)
                arr(i,:)=arr(j,:)
                arr(j,:)=tmp
            end if
            !!按y值
            !if(arr(j,2)<arr(i,2))then
            !    tmp=arr(i,:)
            !    arr(i,:)=arr(j,:)
            !    arr(j,:)=tmp
            !end if
        end do
    end do
end subroutine Bubble_Sort