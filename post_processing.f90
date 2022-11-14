subroutine post_processing
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i

    !创建文件夹
    call makedir

    !打印数值解
    call print_num_data(T_temp)
    
    !打印准确解
    if(case_comp == equEntropy_case .or. case_comp == SodShockTube_case .or. case_comp == LaxShockTube_case .or. case_comp == SinWave_2D_case)then
        call print_exa_data(T_temp)
        !计算范数误差
        call Norm
    end if

    !打印本次计算条件和一些结果
    call print_result
    
    !清除数组内存                 
    do i = 1,ncells+nbdsides
        deallocate( cellset(i).sp_coor,cellset(i).MJacobi,cellset(i).Mdirect,cellset(i).det_J,cellset(i).spvalue_ori,&
                    cellset(i).spvalue_con,cellset(i).spvalue_fluF,cellset(i).spvalue_fluG,&
                    cellset(i).spvalue_con_tem,cellset(i).spvalue_ori_exa,cellset(i).Beta_line)
        deallocate( cellset(i).spvalue_con_loc,cellset(i).spvalue_fluF_loc,cellset(i).spvalue_fluG_loc)      
        deallocate( cellset(i).fluxF_innerfp,cellset(i).fluxG_innerfp)
        deallocate( cellset(i).fpMdirect_F,cellset(i).fpMdirect_G,cellset(i).fpdet_J_F,cellset(i).fpdet_J_G)
        deallocate( cellset(i).spvalue_ori_old)
        deallocate( cellset(i).Smooth_line_x,cellset(i).Smooth_line_y)
    end do

    do i = 1,nsides
        deallocate( sideset(i).fpvalue_upw)
    end do
    deallocate(xy_coor,cellset,sideset,SPs_local,GPs,GCoe,SPs,LPs,FPs,gl_coe,gr_coe,dis_sp_fp,BoundCells_index_set,K_La_Le,K_La_Le_2)
    if(grid_set == self)then
        deallocate(xx,yy)
    end if
end subroutine post_processing

subroutine makedir
    !Ref: Fortran代码自动创建文件夹 https://blog.csdn.net/beyond_buaa/article/details/104168940
    use IFPORT
    implicit none
    integer(kind=4) :: istatus1,errnum
    logical(kind=4) :: ierr1
    character(len=100) :: path_create,path_create2,path_create3
    
    path_create = '.\Result'
    path_create2 = '.\Result\Exact'
    path_create3 = '.\error'
    
    inquire(DIRECTORY=trim(adjustl(path_create)), EXIST=ierr1)
    !print*,ierr1
    if(ierr1) then
       !print*,'The directory have existed and not been needed to create'
       write(*,'(/)')
    else
       istatus1 = SYSTEM('md '//trim(adjustl(path_create)))
       if(istatus1 == -1) then
           errnum=ierrno()
           print*,'Error=',errnum,' inquire the Intel Visual Fortran help document'
           stop ' Folder creating is fail'
       end if
    end if
    
    inquire(DIRECTORY=trim(adjustl(path_create2)), EXIST=ierr1)
    if(ierr1) then

    else      
       istatus1 = SYSTEM('md '//trim(adjustl(path_create2)))
    end if    
    
    inquire(DIRECTORY=trim(adjustl(path_create3)), EXIST=ierr1)
    if(ierr1) then
       
    else      
       istatus1 = SYSTEM('md '//trim(adjustl(path_create3)))
    end if 
end subroutine makedir
subroutine print_submesh_data
    use global_var
    use parameter_setting
    use type_module
    implicit none
    integer :: i,j,k,l,m,P1,P2,P3,P4,iL,iR,LC_sideth,RC_sideth,rk
    character(len=100)filename
    character(len=5)char_nsdpx,char_nsdpy,char_nsdx,char_nsdy
    integer,allocatable ::  SC_BP_index(:,:) !SubCell_Bound_Points_Index(4,nsp),取出靠近单元侧边的点
    integer :: Vertex_P_index(4)
    integer :: cells_contain_vertex(10,2),sum_cells_contain_vertex,nextNearSideth,indexNearCell
    integer :: count1,count2,count3,count4,startCell,nextCell,vertex_th
    real(prec) :: d_D,d_R,d_U,d_L

    open(21,status='REPLACE',file = 'temp.dat')!子单元这些编号可以先存入临时文件。待正式文件写入子单元数后，由临时文件读到正式文件里
    !单元内部解点构成四边形，解点在单元内部编号，从下到上，从左到右:   1,         2,         ..., (j-1)*nsp
    !                                                                  ...
    !                                                                 (nsp-1)*(nsp)+1, j*(nsp)+2, ...,  nsp*nsp
    !write(21,*) 'ZONE, N=',ncells*nsp*nsp,', E=',ncells*(nsp-1)*(nsp-1)+(nsides-nbdsides)*(nsp-1),', F=FEPOINT, ET=QUADRILATERAL'
    !do i = 1,ncells
    !    do j = 1,nsp
    !        do k = 1,nsp
    !            write(21,'(2F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2)
    !        end do
    !    end do
    !end do     
    cycle1:do i = 1,ncells
        cycle2:do j = 1,nsp-1
            !if(j==nsp)cycle
            cycle3:do k = 1,nsp-1  
                P1 = (i-1)*nsp**2+(j-1)*nsp+k
                P2 = P1+1
                P3 = P2+nsp
                P4 = P3-1
                write(21,*) P1,P2,P3,P4
            end do cycle3
        end do cycle2
    end do cycle1

    !do i = 1,nsides
    !    write(*,*)sideset(i).nearcells
    !end do

    !以下是根据单元侧边，连接相对的解点，构成四边形
    allocate(SC_BP_index(4,nsp))
    do j = 1,nsp!从左到右，从下到上
        SC_BP_index(1,j) = j
        SC_BP_index(2,j) = j*nsp
        SC_BP_index(3,j) = (nsp-1)*nsp+j
        
        SC_BP_index(4,j) = (j-1)*nsp+1
    end do

    !以下是靠近单元顶点的解点编号
    Vertex_P_index(1) = 1
    Vertex_P_index(2) = nsp
    Vertex_P_index(3) = nsp*nsp
    Vertex_P_index(4) = (nsp-1)*nsp+1
    do i = 1,nsides
        iL = sideset(i).nearcells(1)
        iR = sideset(i).nearcells(2)    
        if(iL==0 .OR. iR == 0)then
            cycle
        else
            do j = 1,4
                if(cellset(iL).nearcells(j)==iR)then
                    LC_sideth = j !Left cell Side -th.此边是该左单元第j条侧边
                end if
                if(cellset(iR).nearcells(j)==iL)then
                    RC_sideth = j
                end if
            end do
            do k = 1,nsp-1
                P1 = SC_BP_index(LC_sideth,k)+(iL-1)*nsp*nsp
                P2 = SC_BP_index(LC_sideth,k+1)+(iL-1)*nsp*nsp
                if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then
                    rk = nsp+1-k
                    P3 = SC_BP_index(RC_sideth,rk-1)+(iR-1)*nsp*nsp
                    P4 = SC_BP_index(RC_sideth,rk)+(iR-1)*nsp*nsp
                else
                    P3 = SC_BP_index(RC_sideth,k+1)+(iR-1)*nsp*nsp
                    P4 = SC_BP_index(RC_sideth,k)+(iR-1)*nsp*nsp
                end if
                write(21,*) P1,P2,P3,P4
            end do     
        end if
    end do
    !以下是单元顶点处的连接。但此处有点复杂，可能存在三角形，四边形，五边形，六边形，甚至更多边形。后处理的时候要补足空白，使流场中间不留白。
    !   根据顶点查找，然后取出顶点周围单元的最近的解点。连接成四边形+三角形的形式
    !write(*,*) ncells,nbdsides  
    
    sum_vertexCell = 0
    Loop1:do i = 1,nnodes    
        !write(*,*)i
        cells_contain_vertex = 0
        Loop2:do j = 1,ncells
            Loop3:do k = 1,4
                if(cellset(j).nodes(k)==i)then!查询出包含此节点的第一个单元，编号为j
                    startCell = j
                    exit Loop2
                end if
            end do Loop3
        end do Loop2

        ! 之后求邻单元。
        ! 顶点是否在边界上。is: count1 = 1; no: count1 =  0
        call if_VertexInBound(i,count1)
        
        if(count1 == 1)then     !count1 = 1，说明此顶点存在于计算域边界上
            !接下来查询出包含此顶点的所有单元，利用已知的单元可以从其想邻单元入手，缩小查询范围  
            !两侧边都需要查询，因为到边界停止，不确定第一个单元所处位置
            !第一侧边
            !write(*,*)i,k,startCell
            nextNearSideth = k-1
            if(nextNearSideth == 0)nextNearSideth=4
            nextCell = cellset(startCell).nearcells(nextNearSideth)
            !write(*,*)cellset(j).nearcells,nextCell
            
            count2 = 1                       
            cells_contain_vertex(count2,1) = startCell!记录第一个单元
            cells_contain_vertex(count2,2) = k        !记录第几个顶点
            do while((nextCell < ncells+1) .and. (nextCell > 0))   
                !write(*,*)nextCell
                do m = 1,4
                    if(cellset(nextCell).nodes(m)==i) vertex_th = m!m为该顶点在nextCell中的编号
                end do
  
                nextNearSideth = vertex_th-1
                if(nextNearSideth == 0)nextNearSideth=4
                if(cellset(nextCell).nearcells(nextNearSideth)==startCell)then!排除startCell这个nextCell的邻单元
                    nextNearSideth = vertex_th
                    startCell = nextCell
                    nextCell = cellset(startCell).nearcells(nextNearSideth)                    
                elseif(cellset(nextCell).nearcells(vertex_th)==startCell)then
                    startCell = nextCell
                    nextCell = cellset(startCell).nearcells(nextNearSideth)
                else
                    nextCell = 0
                end if                          
                count2 = count2+1
                cells_contain_vertex(count2,1) = startCell!记录单元
                cells_contain_vertex(count2,2) = vertex_th
            end do
            !第二侧边
            nextNearSideth = k
            startCell = j
            nextCell = cellset(startCell).nearcells(nextNearSideth)
            do while((nextCell < ncells+1) .and. (nextCell > 0))   
                !write(*,*)nextCell
                
                do m = 1,4
                    if(cellset(nextCell).nodes(m)==i)vertex_th = m!m为该顶点在nextCell中的编号
                end do
                nextNearSideth = vertex_th-1
                if(nextNearSideth == 0)nextNearSideth=4
                if(cellset(nextCell).nearcells(nextNearSideth)==startCell)then!排除startCell这个nextCell的邻单元
                    nextNearSideth = vertex_th
                    startCell = nextCell
                    nextCell = cellset(startCell ).nearcells(nextNearSideth)
                elseif(cellset(nextCell).nearcells(vertex_th)==startCell)then
                    startCell = nextCell
                    nextCell = cellset(startCell).nearcells(nextNearSideth)
                else
                    nextCell = 0
                end if                        
                count2 = count2+1
                cells_contain_vertex(count2,1) = startCell!记录单元
                cells_contain_vertex(count2,2) = vertex_th
            end do
        elseif(count1 == 0)then     !计算域内部的顶点，可能包含于三个及以上单元中
            !k是该顶点在第一个单元的编号。由此可以知道其k-1,k条侧边所对的相邻单元也包含此顶点。k=1时，k-1设为4
            !利用while循环一直查询邻单元，直到指向第一个单元       
            nextNearSideth = k-1
            if(nextNearSideth == 0)nextNearSideth=4
            nextCell = cellset(startCell).nearcells(nextNearSideth)
            count2 = 1                      
            cells_contain_vertex(count2,1) = startCell!记录第一个单元
            cells_contain_vertex(count2,2) = k
            do while(nextCell /= j)                           

                do m = 1,4
                    if(cellset(nextCell).nodes(m)==i)vertex_th = m!m为该顶点在nextCell中的编号
                end do
 
                nextNearSideth = vertex_th-1
                if(nextNearSideth == 0)nextNearSideth=4
                if(cellset(nextCell).nearcells(nextNearSideth)==startCell)then!排除startCell这个nextCell的邻单元
                    nextNearSideth = vertex_th
                    startCell = nextCell
                    nextCell = cellset(startCell).nearcells(nextNearSideth)
                elseif(cellset(nextCell).nearcells(vertex_th)==startCell)then
                    startCell = nextCell
                    nextCell = cellset(startCell).nearcells(nextNearSideth)
                end if
                count2 = count2+1
                cells_contain_vertex(count2,1) = startCell
                cells_contain_vertex(count2,2) = vertex_th
            end do
        end if                            
        !write(*,*)i,'-----------------------------------------------'
        !write(*,*)cells_contain_vertex(1:5,1)
        
        count3 = 0
        do j = 1,10
            if(cells_contain_vertex(j,1)/=0)then
                count3 = count3+1   !统计包含该顶点单元数
            end if
        end do
        !write(*,*)'2'
        if(count3 == 3)then  !count3 <= 2是边界上的单元且只有两个单元，不用画出四边形
            P1 = (cells_contain_vertex(1,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1,2))
            P2 = (cells_contain_vertex(2,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(2,2))
            P3 = (cells_contain_vertex(3,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(3,2))
            P4 = (cells_contain_vertex(1,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1,2))
            write(21,*) P1,P2,P3,P4
            sum_vertexCell = sum_vertexCell + 1
        elseif(count3 == 4)then !四个单元
            P1 = (cells_contain_vertex(1,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1,2))
            P2 = (cells_contain_vertex(2,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(2,2))
            P3 = (cells_contain_vertex(3,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(3,2))
            P4 = (cells_contain_vertex(4,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(4,2))
            write(21,*) P1,P2,P3,P4
            sum_vertexCell = sum_vertexCell + 1
        elseif(count3 > 4)then                     
            !4个以上单元.根据排列情况:  p = count3 |  四边形    三角形
            !                           p为奇数    |  p/2 - 1      1                               
            !                           p为偶数    | (p-1)/2 - 1   0    
            !偶数情况和奇数情况可以写成一样的形式
            count4 = (count3+1)/2 - 1   !会取整
            !write(*,*)count3,count4
            do j = 1,count4
                !write(*,*)j
                P1 = (cells_contain_vertex(1+(j-1),1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1+(j-1),2))
                P2 = (cells_contain_vertex(2+(j-1),1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(2+(j-1),2))
                P3 = (cells_contain_vertex(count3-j,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(count3-j,2))
                P4 = (cells_contain_vertex(count3-(j-1),1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(count3-(j-1),2))
                write(21,*) P1,P2,P3,P4
                sum_vertexCell = sum_vertexCell + 1
            end do     
        end if
    end do Loop1

    !write(*,*)sum_vertexCell
    !do i =1,ncells
    !    write(*,*)i,cellset(i).nodes
    !end do
    
    !write(*,*)((SC_BP_index(i,j),j=1,nsp),i=1,4)
    close(21)
    deallocate(SC_BP_index)
end subroutine print_submesh_data
subroutine if_VertexInBound(indexVertex,res)
    use type_module
    use real_precision
    use global_var
    use bc_module
    implicit none
    integer :: indexVertex,indexSide,indexCell,res,i,j
    
    res = 0
    loop1:do i = 1,nbdsides
        indexCell = BoundCells_index_set(i)
        loop2:do j = 1,4
            indexSide = cellset(indexCell).sides(j)
            if(sideset(indexSide).bc /= Interior)then
                if(indexVertex == sideset(indexSide).nodes(1) .or. indexVertex == sideset(indexSide).nodes(2))then
                    res = 1
                    exit loop1
                end if
            end if
        end do loop2
    end do loop1
    
end subroutine if_VertexInBound
    
subroutine print_mesh_data
    use global_var
    use parameter_setting
    use type_module
    use bc_module
    implicit none
    integer :: i,j,k,l,m,P1,P2,P3,P4,iL,iR,LC_sideth,RC_sideth,rk
    character(len=60)filename
    character(len=5)char_nsdpx,char_nsdpy,char_nsdx,char_nsdy
    character(len=20)lines(4)
    integer :: status1
    integer :: Vertex_P_index(4)
    integer :: cells_contain_vertex(10,2),sum_cells_contain_vertex,nextNearSideth,indexNearCell
    integer :: count1,count2,count3,count4,startCell,nextCell,vertex_th
    real(prec) :: d_D,d_R,d_U,d_L
    !----输出.plt查看 Mesh------------------------------------------------------------------------------   
    select case(grid_set)
    case(self)
        write(char_nsdpx,"(TL1,I4)") nsdpx-1
        write(char_nsdpy,"(TL1,I4)") nsdpy-1
        select case(cell_shape)
        case(cell_rect)
            filename = '.\mesh\mesh_rect_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'.plt'
        case(cell_trap)
            filename = '.\mesh\mesh_trap_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'.plt'
        case(cell_dist)
            filename = '.\mesh\mesh_dist_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'.plt'
        end select   
    case(fluent_cas)
        write(char_nsdx,"(TL1,I4)") nsdx
        write(char_nsdy,"(TL1,I4)") nsdy
        filename = '.\mesh\mesh_'//trim(mesh_file)//'.plt'
        
    end select
    
    open(20,status='REPLACE',file = filename)
    write(20,*)'title="Mesh"'
    write(20,*) 'variables =x,y'
    write(20,*) 'ZONE, N=',nnodes,', E=',ncells,', F=FEPOINT, ET=QUADRILATERAL'
    
    !FEPOINT形式
    do i = 1,nnodes
        write(20,*) xy_coor(i,:)
        !write(*,*) xy_coor(i,:)
    end do
    !四边形单元顶点编号
    do i = 1,ncells
        write(20,*) cellset(i).nodes(1:4)
    end do
       
    write(20,*) 'ZONE, N=',ncells*nsp*nsp,', E=',ncells*(nsp-1)*(nsp-1)+(nsides-nbdsides)*(nsp-1)+sum_vertexCell,', F=FEPOINT, ET=QUADRILATERAL'   
    !解点子单元值
    do i = 1,ncells
        do j = 1,nsp
            do k = 1,nsp
                write(20,'(2F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2)
            end do
        end do
    end do 
    !解点子单元编号
    open(UNIT = 21 , FILE = 'temp.dat' , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1)  ! 打开文件
    do while(status1 == 0 )
        read(21 , * , IOSTAT = status1) lines               
        write(20,'(4A20)')lines
    end do
    
    close(21)
    close(20)
end subroutine print_mesh_data

subroutine print_num_data(T_cal)
    use global_var
    use parameter_setting
    use type_module
    use preParameter
    implicit none
    integer :: i,j,k,l,P1,P2,P3,P4,iL,iR,LC_sideth,RC_sideth,rk
    character(len=100)filename
    character(len=16)char_nsdpx,char_nsdpy,char_T,char_nsp,char_case,char_dire,char_detect,char_buffer,char_Dim,char_scheme_kind,char_Ms,char_Mv
    character(len=20)lines(4)
    integer :: status1,counter_beta
    real(prec) :: percentage,T_cal

    !----输出.plt查看 Result------------------------------------------------------------------------------   
    !数值解

    !创建文件夹
    call makedir
    
    select case(detection_type)
    case(detection_TVB)
        char_detect = 'TVB'
    case(detection_ATV)
        char_detect = 'ATV'
    case(detection_MV)       
        char_detect = 'MV'
    case(detection_MDH)
        char_detect = 'MDH'
    case(detection_MDHm)
        char_detect = 'MDHm'
    case(detection_KXRCF)
        char_detect = 'KXRCF'
    case(detection_JST)
        char_detect = 'JST'
    case default
        char_detect = 'newadd'
    end select
    
    select case(detect_type)
    case(ByDim )
        char_Dim = 'ByDim'
    case(ByCell)
        char_Dim = 'ByCell'
    end select
    
    select case(scheme_kind)
    case(scheme_cpr)
        char_scheme_kind = 'CPR'
    case(scheme_two)
        char_scheme_kind = 'Two'
    case(scheme_hybrid)
        char_scheme_kind = 'hybrid'
    end select
    
    write(char_T,"(TL1,F16.4)") T
    if(T_cal<T)then
        write(char_T,"(TL1,F16.4)") T_cal
    end if
    write(char_case,"(TL1,I2)")case_comp
    
    select case(grid_set)
    case(self)
        write(char_nsdpx,"(TL1,I4)") nsdpx-1
        write(char_nsdpy,"(TL1,I4)") nsdpy-1
        write(char_nsp,"(TL1,I4)")nsp
        
        write(char_dire,"(TL1,I1)")dire_shock

        select case(cell_shape)
        case(cell_rect)
            filename = '.\Result\'//'Z_case'//trim(adjustl(char_case))//'__'//trim(adjustl(char_nsdpx))//'_'//trim(adjustl(char_nsdpy))//'_nsp_'//&
            trim(adjustl(char_nsp))//'_T_'//trim(adjustl(char_T))//'_'//trim(adjustl(char_detect))//trim(adjustl(char_Dim))//'_'//trim(adjustl(char_scheme_kind))
        case(cell_trap)
            filename = '.\Result\'//'Zhe_case'//trim(adjustl(char_case))//'__'//trim(adjustl(char_nsdpx))//'_'//trim(adjustl(char_nsdpy))//'_nsp_'//&
            trim(adjustl(char_nsp))//'_T_'//trim(adjustl(char_T))//'_'//trim(adjustl(char_detect))//trim(adjustl(char_Dim))//'_'//trim(adjustl(char_scheme_kind))
        case(cell_dist)
            filename = '.\Result\'//'Dis_case'//trim(adjustl(char_case))//'__'//trim(adjustl(char_nsdpx))//'_'//trim(adjustl(char_nsdpy))//'_nsp_'//&
            trim(adjustl(char_nsp))//'_T_'//trim(adjustl(char_T))//'_'//trim(adjustl(char_detect))//trim(adjustl(char_Dim))//'_'//trim(adjustl(char_scheme_kind))
        end select  
        if(case_comp == CompositeVortexShock_case)then
            write(char_Ms,"(TL1,F7.2)") Ms
            write(char_Mv,"(TL1,F7.2)") Mv
            filename = trim(adjustl(filename))//'_Mv'//trim(adjustl(char_Mv))//'_Ms'//trim(adjustl(char_Ms))
        end if
        if(buffer_cell_switch==buffer_cell_yes)then
            filename = trim(adjustl(filename))//'_buffer'
        end if        
        filename = trim(adjustl(filename)) //'.plt'
        
        open(40,status= 'REPLACE',file = filename)
        
        write(40,*)'variables =x,y,r,u,v,p,Beta,BetaX,BetaY,BetaXY,alpha_x,alpha_y'
        !结构形式
        write(40,*)'zone i=',nsp*(nsdpx-1),'j=',nsp*(nsdpy-1),'   F = point' 
        !write(40,*)'zone i=',nsp*(nsdpx-1),'j=',nsp*(nsdpy-1)+nsp,'   F = point'!输出最上层时需要
        do l = 1,nsdy
            do j = 1,nsp 
                do i = (l-1)*(nsdpx-1)+1,l*(nsdpx-1)
                    do k =1,nsp            
                        write(40,'(6F22.17,4I4,2F22.17)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori(j,k,1),cellset(i).spvalue_ori(j,k,2),&
                                                         cellset(i).spvalue_ori(j,k,3),cellset(i).spvalue_ori(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2),&
                                                         cellset(i).Beta_line(j,1)+cellset(i).Beta_line(k,2),cellset(i).Smooth_line_x(j,k),cellset(i).Smooth_line_y(k,j)
                        !write(40,'(6F22.17I4I4I4)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_con(j,k,1),cellset(i).spvalue_con(j,k,2),cellset(i).spvalue_con(j,k,3),cellset(i).spvalue_con(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2)
                    end do
                end do   
            end do
        end do

        !do j = 1,nsp 
        !    do i = ncells+n_BCells_D+ n_BCells_R+1,ncells+n_BCells_D+ n_BCells_R + n_BCells_U
        !        do k =1,nsp            
        !            write(40,'(6F22.17,4I4)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori(j,k,1),cellset(i).spvalue_ori(j,k,2),cellset(i).spvalue_ori(j,k,3),cellset(i).spvalue_ori(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2),cellset(i).Beta_line(j,1)+cellset(i).Beta_line(k,2)
        !            !write(40,'(6F22.17I4I4I4)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_con(j,k,1),cellset(i).spvalue_con(j,k,2),cellset(i).spvalue_con(j,k,3),cellset(i).spvalue_con(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2)
        !        end do
        !    end do   
        !end do
        
        !!非结构形式
        !write(40,*) 'ZONE, N=',ncells*nsp*nsp,', E=',ncells*(nsp-1)*(nsp-1)+(nsides-nbdsides)*(nsp-1)+sum_vertexCell,', F=FEPOINT, ET=QUADRILATERAL' 
        !解点子单元值
        ! do i = 1,ncells
        !    do j = 1,nsp
        !        do k = 1,nsp
        !            write(40,'(6F20.10,4I4)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori(j,k,1),cellset(i).spvalue_ori(j,k,2),cellset(i).spvalue_ori(j,k,3),cellset(i).spvalue_ori(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2),cellset(i).Beta_line(j,1)+cellset(i).Beta_line(k,2)
        !        end do
        !    end do
        !end do    
        !!解点子单元编号
        !open(UNIT = 21 , FILE = 'temp.dat' , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1)  ! 打开文件
        !do while(status1 == 0 )
        !    read(21 , * , IOSTAT = status1) lines               
        !    write(40,'(4A20)')lines
        !end do 

    case(Fluent_cas)

        filename = '.\Result\Unstru_case'//trim(adjustl(char_case))//'__T'//trim(adjustl(char_T))//'_'//trim(adjustl(char_detect))//trim(adjustl(char_Dim))//'_'//trim(mesh_file)//'_'//trim(adjustl(char_scheme_kind))//'.plt'
        open(40,status='REPLACE',file = filename)
                
        write(40,*)'variables =x,y,r,u,v,p,Beta,BetaX,BetaY,BetaXY,alpha_x,alpha_y'
        write(40,*) 'ZONE, N=',ncells*nsp*nsp,', E=',ncells*(nsp-1)*(nsp-1)+(nsides-nbdsides)*(nsp-1)+sum_vertexCell,', F=FEPOINT, ET=QUADRILATERAL'   
        !解点子单元值
         do i = 1,ncells
            do j = 1,nsp
                do k = 1,nsp
                    write(40,'(6F22.17,4I4,2F22.17)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori(j,k,1),cellset(i).spvalue_ori(j,k,2),&
                                                     cellset(i).spvalue_ori(j,k,3),cellset(i).spvalue_ori(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2),&
                                                     cellset(i).Beta_line(j,1)+cellset(i).Beta_line(k,2),cellset(i).Smooth_line_x(j,k),cellset(i).Smooth_line_y(k,j)
                end do
            end do
        end do
        !解点子单元编号
        open(UNIT = 21 , FILE = 'temp.dat' , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1)  ! 打开文件
        do while(status1 == 0 )
            read(21 , * , IOSTAT = status1) lines               
            write(40,'(4A20)')lines
        end do    
    end select
    
    
    !打印壁面信息.只是绕流   
    if(case_comp == WingFlow_case .or. case_comp == HyperCylinder_case )then
        call print_wall_data
    end if
      
    !输入计算参数
    write(40,NML=PARA)
    write(40,NML=PARA_debug) 
    counter_beta = 0
    do i = 1, ncells
        if(cellset(i).Beta==1)then
            counter_beta = counter_beta + 1
        end if
    end do
    
    percentage = real(counter_beta)/ncells
    write(*,"('TCs : 'F10.6)")percentage
    write(40,"('Trouble cells percentage: 'F10.6)")percentage
    
    close(21)
    
    if(compuPurpose == error_solve)then
        
    else
        close(40)
    end if
    !close(40)
end subroutine print_num_data

subroutine print_exa_data(T_cal)
    use global_var
    use parameter_setting
    use type_module
    implicit none
    integer :: i,j,k,l,P1,P2,P3,P4
    character(len=100)filename
    character(len=10)char_nsdpx,char_nsdpy,char_T,char_nsp,char_case,char_dire
    character(len=20)lines(4)
    integer :: status1
    real(prec) :: T_cal
    !----输出.plt查看 Result------------------------------------------------------------------------------   
    !准确解
    select case(grid_set)
    case(self)
        write(char_nsdpx,"(TL1,I4)") nsdpx-1
        write(char_nsdpy,"(TL1,I4)") nsdpy-1
        write(char_nsp,"(TL1,I4)")nsp
        write(char_case,"(TL1,I4)")case_comp
        write(char_dire,"(TL1,I1)")dire_shock
        write(char_T,"(TL1,F8.4)") T_cal
        if(nt_temp*dt>T)then
            write(char_T,"(TL1,F6.4)") T_cal
        end if
        select case(cell_shape)
        case(cell_rect)
            filename = '.\Result\Exact\Z_case_'//trim(adjustl(char_case))//'_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'_nsp_'//&
            trim(adjustl(char_nsp))//'_T_'//trim(adjustl(char_T))//'_'//trim(adjustl(char_dire))//'.plt'
        case(cell_trap)
            filename = '.\Result\Exact\Zhe_case_'//trim(adjustl(char_case))//'_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'_nsp_'//&
            trim(adjustl(char_nsp))//'_T_'//trim(adjustl(char_T))//'_'//trim(adjustl(char_dire))//'.plt'
        case(cell_dist)
            filename = '.\Result\Exact\Dis_case_'//trim(adjustl(char_case))//'_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'_nsp_'//&
            trim(adjustl(char_nsp))//'_T_'//trim(adjustl(char_T))//'_'//trim(adjustl(char_dire))//'.plt'
        end select 
        
        open(60,file = filename)
        write(60,*)'variables =x,y,r,u,v,p'
        write(60,*)'zone i=',nsp*(nsdpx-1),'j=',nsp*(nsdpy-1),'   F = point' 
    
        do l = 1,nsdy
            do j = 1,nsp 
                do i = (l-1)*(nsdpx-1)+1,l*(nsdpx-1) 
                    do k =1,nsp  
                        !write(*,'(6F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori_exa(j,k,1),cellset(i).spvalue_ori_exa(j,k,2),cellset(i).spvalue_ori_exa(j,k,3),cellset(i).spvalue_ori_exa(j,k,4)               
                        write(60,'(6F22.17)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori_exa(j,k,1),cellset(i).spvalue_ori_exa(j,k,2),cellset(i).spvalue_ori_exa(j,k,3),cellset(i).spvalue_ori_exa(j,k,4)
                        !write(60,'(6F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_con(j,k,1),cellset(i).spvalue_con(j,k,2),cellset(i).spvalue_con(j,k,3),cellset(i).spvalue_con(j,k,4)
                    end do
                end do   
            end do
        end do
    case(fluent_cas)
        write(char_T,"(TL1,F6.4)") T
        write(char_case,"(TL1,I1)")case_comp
        filename = '.\Result\Exact\Unstru_case'//trim(adjustl(char_case))//'_'//trim(mesh_file)//'_T_'//trim(adjustl(char_T))//'.plt'
        open(60,file = filename)
        write(60,*)'variables =x,y,r,u,v,p'    
        write(60,*) 'ZONE, N=',ncells*nsp*nsp,', E=',ncells*(nsp-1)*(nsp-1)+(nsides-nbdsides)*(nsp-1)+sum_vertexCell,', F=FEPOINT, ET=QUADRILATERAL'   
        !解点子单元值
         do i = 1,ncells
            do j = 1,nsp
                do k = 1,nsp
                    write(60,'(6F22.17)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori_exa(j,k,1),cellset(i).spvalue_ori_exa(j,k,2),cellset(i).spvalue_ori_exa(j,k,3),cellset(i).spvalue_ori_exa(j,k,4)
                end do
            end do
        end do    
        !解点子单元编号
        open(UNIT = 21 , FILE = 'temp.dat' , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1)  ! 打开文件
        do while(status1 == 0 )
            read(21 , * , IOSTAT = status1) lines               
            write(60,'(4A20)')lines
        end do    
    end select
    
    close(60)
    close(21)
end subroutine print_exa_data

subroutine print_result
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,counter_beta
    character(len=100) :: filename
    character(len=10) :: char_detect
    real(prec) :: percentage
    filename = 'Result of this calculation.txt'
    open(1100,status='REPLACE',file = filename)
    write(1100,NML=PARA)
    write(1100,NML=PARA_debug)
    counter_beta = 0
    do i = 1, ncells
        if(cellset(i).Beta==1)then
            counter_beta = counter_beta + 1
        end if
    end do
    
    select case(detection_type)
    case(detection_TVB)
        char_detect = 'TVB'
    case(detection_ATV)
        char_detect = 'ATV'
    case(detection_MV)       
        char_detect = 'MV'
    case(detection_MDH)
        char_detect = 'MDH'
    case(detection_MDHm)
        char_detect = 'MDHm'
    case(detection_KXRCF)
        char_detect = 'KXRCF'
    case(detection_JST)
        char_detect = 'JST'
    case default
        char_detect = 'newadd'
    end select
    
    write(*,"('TCs :',F10.6)")percentage
    percentage = real(counter_beta)/ncells
    write(1100,"('computation time: 'F10.6)")nt_temp*dt
    write(1100,"('trouble cells indicator: 'A10)")char_detect
    write(1100,"('trouble cells percentage: 'F10.6)")percentage
    
    close(1100)
end subroutine print_result

subroutine print_wall_data    
    use global_var
    use parameter_setting
    use type_module
    use bc_module
    implicit none
    integer :: i,j,k,l,m,P1,P2,P3,P4,iL,iR,LC_sideth,RC_sideth,rk
    character(len=100)filename
    character(len=5)char_nsdpx,char_nsdpy,char_nsdx,char_nsdy
    integer,allocatable ::  SC_BP_index(:,:) !SubCell_Bound_Points_Index(4,nsp),取出靠近单元侧边的点
    integer :: Vertex_P_index(4),count_wall_sides
    integer :: cells_contain_vertex(10,2),sum_cells_contain_vertex,nextNearSideth,indexNearCell,sideIndex,cellIndex
    integer :: count1,count2,count3,count4,startCell,nextCell,vertex_th,count_temp
    real(prec),dimension(:,:) :: bc_nodes_values(nsp,4)
    real(prec),dimension(:,:),allocatable :: wallSPs_set
    real(prec) :: tmp(6)
    write(40,*) '' 
    write(40,*) 'zone' 
    !简单寻找出边界附近的求解点
    count_wall_sides = 0
    do i = 1,nbdsides
        cellIndex = BoundCells_index_set(i)
        if(i>1 .and. cellIndex == BoundCells_index_set(i-1))cycle !有拐角单元
        do j = 1,4
            sideIndex = cellset(cellIndex).sides(j)
            if(sideset(sideIndex).bc==Wall)then
                count_wall_sides = count_wall_sides + 1
            end if
        end do
    end do
    allocate(wallSPs_set(count_wall_sides*nsp,6))
 
    count_temp = 1
    do i = 1,nbdsides
        cellIndex = BoundCells_index_set(i)
        if(i>1 .and. cellIndex == BoundCells_index_set(i-1))cycle !有拐角单元
        do j = 1,4
            sideIndex = cellset(cellIndex).sides(j)
            if(sideset(sideIndex).bc==Wall)then
                if(j == 1)then                     
                    do k = 1,nsp
                        wallSPs_set(count_temp,1:2) = cellset(cellIndex).sp_coor(1,k,:)
                        wallSPs_set(count_temp,3:6) = cellset(cellIndex).spvalue_ori(1,k,:)
                        count_temp  = count_temp + 1
                    end do  
                elseif(j == 2)then
                    do k = 1,nsp
                        wallSPs_set(count_temp,1:2) = cellset(cellIndex).sp_coor(k,nsp,:)
                        wallSPs_set(count_temp,3:6) = cellset(cellIndex).spvalue_ori(k,nsp,:)
                        count_temp  = count_temp + 1
                    end do
                    !write(*,*)cellset(cellIndex).sp_coor(1,nsp,:)
                elseif(j == 3)then
                    do k = 1,nsp
                        wallSPs_set(count_temp,1:2) = cellset(cellIndex).sp_coor(nsp,k,:)
                        wallSPs_set(count_temp,3:6) = cellset(cellIndex).spvalue_ori(nsp,k,:)
                        count_temp  = count_temp + 1
                    end do
                    !write(*,*)cellset(cellIndex).sp_coor(nsp,1,:)
                elseif(j == 4)then
                    do k = 1,nsp
                        wallSPs_set(count_temp,1:2) = cellset(cellIndex).sp_coor(k,1,:)
                        wallSPs_set(count_temp,3:6) = cellset(cellIndex).spvalue_ori(k,1,:)
                        count_temp  = count_temp + 1
                    end do
                    !write(*,*)cellset(cellIndex).sp_coor(1,1,:)
                end if
            end if
        end do
    end do
    
    !排序    
    do i = 1, count_wall_sides*nsp
        do j = i+1,count_wall_sides*nsp
            !!按x值 y>0时
            if(wallSPs_set(j,1) < wallSPs_set(i,1) .and. wallSPs_set(j,2) .GE. 0 .and.wallSPs_set(i,2) .GE. 0 )then
                tmp = wallSPs_set(i,:)
                wallSPs_set(i,:) = wallSPs_set(j,:)
                wallSPs_set(j,:) = tmp
            end if
        end do
    end do
    do i = 1, count_wall_sides*nsp
        do j = i+1,count_wall_sides*nsp
            !!按x值 y<0时
            if(wallSPs_set(j,1) > wallSPs_set(i,1) .and. wallSPs_set(j,2) < 0 .and. wallSPs_set(i,2) < 0)then
                tmp = wallSPs_set(i,:)
                wallSPs_set(i,:) = wallSPs_set(j,:)
                wallSPs_set(j,:) = tmp
            end if
        end do
    end do
    !打印
    do i = 1, count_wall_sides*nsp
        !write(*,*)wallSPs_set(i,1:2)
        write(40,'(6F22.17,6I4)')wallSPs_set(i,1:6),0,0,0,0,0,0
    end do
    
    !
    !count_temp = 1
    !do i = 1,nbdsides
    !    cellIndex = BoundCells_index_set(i)
    !    do j = 1,4
    !        sideIndex = cellset(cellIndex).sides(j)
    !        if(sideset(sideIndex).bc==Wall)then
    !            if(j == 1)then                    
    !                do k = 1,nsp
    !                write(40,'(6F22.17,6I4)')cellset(cellIndex).sp_coor(1,k,1),cellset(cellIndex).sp_coor(1,k,2),cellset(cellIndex).spvalue_ori(1,k,1),cellset(cellIndex).spvalue_ori(1,k,2),cellset(cellIndex).spvalue_ori(1,k,3),cellset(cellIndex).spvalue_ori(1,k,4),0,0,0,0,0,0
    !                end do  
    !                !write(*,*)cellset(cellIndex).sp_coor(1,1,:)
    !            elseif(j == 2)then
    !                do k = 1,nsp
    !                write(40,'(6F22.17,6I4)')cellset(cellIndex).sp_coor(k,nsp,1),cellset(cellIndex).sp_coor(k,nsp,2),cellset(cellIndex).spvalue_ori(k,nsp,1),cellset(cellIndex).spvalue_ori(k,nsp,2),cellset(cellIndex).spvalue_ori(k,nsp,3),cellset(cellIndex).spvalue_ori(k,nsp,4),0,0,0,0,0,0
    !                end do
    !                !write(*,*)cellset(cellIndex).sp_coor(1,nsp,:)
    !            elseif(j == 3)then
    !                do k = 1,nsp
    !                write(40,'(6F22.17,6I4)')cellset(cellIndex).sp_coor(nsp,nsp+1-k,1),cellset(cellIndex).sp_coor(nsp,nsp+1-k,2),cellset(cellIndex).spvalue_ori(nsp,nsp+1-k,1),cellset(cellIndex).spvalue_ori(nsp,nsp+1-k,2),cellset(cellIndex).spvalue_ori(nsp,nsp+1-k,3),cellset(cellIndex).spvalue_ori(nsp,nsp+1-k,4),0,0,0,0,0,0
    !                end do
    !                !write(*,*)cellset(cellIndex).sp_coor(nsp,1,:)
    !            elseif(j == 4)then
    !                do k = 1,nsp
    !                write(40,'(6F22.17,6I4)')cellset(cellIndex).sp_coor(k,1,1),cellset(cellIndex).sp_coor(k,1,2),cellset(cellIndex).spvalue_ori(k,1,1),cellset(cellIndex).spvalue_ori(k,1,2),cellset(cellIndex).spvalue_ori(k,1,3),cellset(cellIndex).spvalue_ori(k,1,4),0,0,0,0,0,0
    !                end do
    !                !write(*,*)cellset(cellIndex).sp_coor(1,1,:)
    !            end if
    !        end if
    !    end do
    !    
    !end do
    !
    !stop
    
end subroutine print_wall_data      
subroutine renewal_program
    !续算程序--数据记录
    !隔步记录结果
    !续算判断指标：时间是否一致（大于，小于，等于）
    !            ：条件是否一致
    !或只有满足某条件才开始续算：预设值完全一致，只有当前时间步小于预设步时开始续算。
    !   ##文件记录格式##
    !   #   前n行记录条件
    !   #   后续记录数据
    use global_var 
    use parameter_setting 
    use real_precision
    use type_module
    
    integer :: i,j,k,l
    character(len=100) :: filename
    character(len=6)   :: char_nsdpx,char_nsdpy,char_T,char_nsp,char_case,char_dire
    character(len=20)  :: lines(4)
    character(len=1024) :: char_namelist
    integer :: status1

    !write(char_namelist,NML=PARA)
    !write(*,*) trim(char_namelist)

    filename = '.\renewal\ccc.con'       
    open(400,status='REPLACE',file = filename)
    
    
    write(400,NML=PARA)
    write(400,NML=PARA_debug)
    write(400,*)'variables =x,y,r,u,v,p,Beta,BetaX,BetaY,BetaXY'
    
    do i = 1,ncells
        do j = 1,nsp
            do k = 1,nsp
                write(400,'(6F22.17,4I4)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori(j,k,1),cellset(i).spvalue_ori(j,k,2),cellset(i).spvalue_ori(j,k,3),cellset(i).spvalue_ori(j,k,4),cellset(i).Beta,cellset(i).Beta_line(j,1),cellset(i).Beta_line(k,2),cellset(i).Beta_line(j,1)+cellset(i).Beta_line(k,2)
            end do
        end do
    end do
        
    close(400)
    
    
end subroutine renewal_program