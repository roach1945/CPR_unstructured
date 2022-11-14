
subroutine Time_com
    use real_precision
    use Time_module
    implicit none
    real(prec) :: Total,t1,t2,t3,t4
    Total = TimeEnd-TimeStart
    t1 = TimeFrame1-TimeStart
    t2 = TimeFrame2-TimeFrame1
    t3 = TimeFrame3-TimeFrame2
    t4 = TimeEnd-TimeFrame3
    write(*,*)
    write(*,"(1xA)")'computional time'
    write(*,"('Total time:'F10.3' s, ',F10.3,' min')")Total,Total/60.0_prec
    write(*,"('pre_processing:'F10.3' s, percentage:'F5.2)") t1,t1/Total
    write(*,"('time_advancing:'F10.3' s, percentage:'F5.2)") t2,t2/Total
    write(*,"('exact_solution:'F10.3' s, percentage:'F5.2)") t3,t3/Total
    write(*,"('post_processing:'F9.3' s, percentage:'F5.2)") t4,t4/Total
    
    write(40,"('Total time:'F10.3' s, ',F10.3,' min')")Total,Total/60.0_prec
    write(40,"('pre_processing:'F10.3' s, percentage:'F5.2)") t1,t1/Total
    write(40,"('time_advancing:'F10.3' s, percentage:'F5.2)") t2,t2/Total
    write(40,"('exact_solution:'F10.3' s, percentage:'F5.2)") t3,t3/Total
    write(40,"('post_processing:'F9.3' s, percentage:'F5.2)") t4,t4/Total
    write(40,"('detection:'F9.3' s, percentage:'F5.2)")DetectionTime,DetectionTime/Total
    write(*,"('detection:'F9.3' s, percentage:'F7.4)")DetectionTime,DetectionTime/Total
    close(40)!!计算结果打印文件
    
end subroutine Time_com 

subroutine read_preParameter
    use preParameter
    implicit none 

    character(len = 80)::msg        ! 打开文件错误时的异常信息
    integer :: nvals = 0            ! 读取数据的行数
    integer :: status1              ! 打开文件和读取数据的状态信息
    real(prec) :: value             ! 需读取文件中的数据
    
    open(UNIT = 10 , FILE = 'input.in' , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1 , IOMSG = msg)  

    
    if(status1 == 0 ) then                                      
        read(10,NML=prePARA) 
    else                                                        
        write(*,"('Error opening file : IOSTAT = ',I6)") status1
        write(*,*) trim(msg)                                   
        stop
    end if
    close(10) 
    end subroutine read_preParameter
    
    
subroutine time_cal

    !-----------------------------------------------------------------------------
    !
    !   方法：估计子程序的时间花销
    !   描述：如测试时间推进中子程序的时间花销。用来测试，目前删掉了。
    !         call cpu_time(T1)
    !         call subroutine_cost
    !         call cpu_time(T2)
    !   作者：gqShi 
    !   状态：未启用
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use time_test
    implicit none
    write(*,*) T1-Ts
    write(*,*) T2-T1
    write(*,*) T3-T2
    write(*,*) T4-T3
    write(*,*) T5-T4
    write(*,*) Te-T5
    pause
    
end subroutine time_cal
    
subroutine correct_density_pressure

    !-----------------------------------------------------------------------------
    !
    !   方法：矫正压力密度出负
    !   描述：每一时间步计算完毕，检测流场中是否出现计算出负的情况，然后根据周边数据进行矫正
    !         取周围4个点的平均值；但效果不好，暂未启用
    !   作者：gqShi 
    !   状态：未启用
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer :: i,j,k,l
    real(prec) :: ori,r,p
    integer :: NearPs(4,3)
    do i = 1,ncells  
        do k = 1,nsp
            do l = 1,nsp
                r = cellset(i).spvalue_ori(k,l,1)
                p = cellset(i).spvalue_ori(k,l,4)                
                if(r< 0.0_prec .or. p< 0.0_prec )then
                    call searchNearPoints(i,k,l,NearPs)
                    cellset(i).spvalue_ori(k,l,:) = (cellset(NearPs(1,1)).spvalue_ori(NearPs(1,2),NearPs(1,3),:)+cellset(NearPs(2,1)).spvalue_ori(NearPs(2,2),NearPs(2,3),:) &
                       +cellset(NearPs(3,1)).spvalue_ori(NearPs(3,2),NearPs(3,3),:)+cellset(NearPs(4,1)).spvalue_ori(NearPs(4,2),NearPs(4,3),:))*0.25_prec
                end if
            end do   
        end do
    end do       

end subroutine correct_density_pressure

subroutine searchNearPoints(index,tk,tl,NearPs)
    use real_precision
    use parameter_setting
    use global_var
    use type_module
    implicit none
    integer :: index,nearIndex,i,j,k,l,tk,tl,NearPs(4,3)
    real(prec) :: temp_dis(4),allNearPs(5,nsp,nsp,2)!(4点;i，k，l)
    real(prec) ::dis
    allNearPs(1,:,:,:) = cellset(index).sp_coor(:,:,:)
    !write(*,*)'000',index    
    do i = 2,5
        NearIndex = cellset(index).nearcells(i-1)
        write(*,*)NearIndex
        allNearPs(i,:,:,:) = cellset(NearIndex).sp_coor(:,:,:)
    end do
    write(*,*)'001'
    temp_dis = 100.0_prec
    do i = 1,5
        write(*,*)i
        do k = 1,nsp
            do l =1,nsp
                if(k == tk .and. l == tl)then
                    cycle
                end if
                dis = (allNearPs(i,k,l,1)-allNearPs(1,tk,tl,1))**2 + (allNearPs(i,k,l,2)-allNearPs(1,tk,tl,2))**2
                !write(*,*)'00',dis
                do j = 1,4
                    if(dis<temp_dis(j))then
                        temp_dis(j)=dis
                        NearPs(j,2) = k
                        NearPs(j,3) = l
                        !write(*,*)'02',j,k,l
                        if(i==1)then
                            NearPs(j,1) = index
                        else
                            NearPs(j,1) = cellset(index).nearcells(i-1)
                        end if 
                        exit
                    end if
                end do
            end do
        end do
    end do
    write(*,*)NearPs
end subroutine searchNearPoints
    
subroutine com_time(n)
    
    !-----------------------------------------------------------------------------
    !
    !   方法：估计时间推进计算时间
    !   描述：计算20步，然后用总步数nstep估计总时间。若是自适应时间步长，预估有可能相差很大，固定时间步长相差不大。
    !         可以更改 21 -> another number 以改变用来评估的步数
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use parameter_setting
    use global_var,only:hour_whole,minite_whole,second_whole
    implicit none
    character(8) :: date
    character(10) :: time
    character(5) :: zone
    integer :: values(8),hh,mm,ss,time_long,n,dh,dm,ds,oh,om,os

    call DATE_AND_TIME(TIME = time)
    read(time,"(I2,I2,I2)") hh,mm,ss
   
    if(n==1)then
        hour_whole = hh
        minite_whole = mm
        second_whole = ss
    end if
  
    if(n == 21)then
        write(*,"('    ### Now Time : 'I2'h 'I2'm 'I2's ###')") hh,mm,ss
        hour_whole = hh - hour_whole
        minite_whole = mm - minite_whole
        second_whole = ss - second_whole
        time_long = hour_whole*3600 + minite_whole*60 + second_whole
        time_long = (nstep/20 -1)* time_long
        dh = time_long/3600
        dm = time_long/60 - dh*60
        ds = time_long - dh*3600 - dm*60       
        write(*,"('    ### Com Time : 'I4'h 'I3'm 'I3's ###')")dh,dm,ds
        os = mod(ss+ds,60)
        om = mod((mm + dm) + (ss + ds)/60,60)
        oh = (hh + dh) + ((mm + dm) + (ss + ds)/60)/60
        write(*,"('    ### Ove Time : 'I4'h 'I3'm 'I3's ###')")oh,om,os
    end if
    
    !write(*,*)hh,mm,ss
end subroutine com_time
    
subroutine now_time

    !-----------------------------------------------------------------------------
    !
    !   方法：打印现在的时间
    !   描述：评估计算进度
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    implicit none
    character(8) :: date
    character(10):: time
    character(5) :: zone
    integer :: values(8),hh,mm,ss

    call DATE_AND_TIME(TIME = time)
    
    read(time,"(I2,I2,I2)") hh,mm,ss
    write(*,"('  Now Time : 'I2'h 'I2'm 'I2's')") hh,mm,ss
    
end subroutine now_time

subroutine solve_errorGCL(ntimestep)

    use real_precision
    use parameter_setting
    use global_var
    use type_module
    implicit none
    
    integer :: i,j,k,l,ntimestep
    real(prec) :: var_array(nsp,nsp),detJ_array(nsp,nsp),sum,sum_temp,detJ_temp,sum_q_exact
    
    sum = 0.0_prec
    
    !数值解的全局密度通量
    do i = 1,ncells
        var_array = 0.0_prec
        sum_temp = 0.0_prec
        var_array = cellset(i).spvalue_ori(:,:,1)*cellset(i).det_J(:,:)
        call Gauss_DInt(var_array,sum_temp)       
        sum = sum + sum_temp
    end do
    sum_q = sum
    if(ntimestep == 0)then
        sum_q_initial = sum
    end if
    
    !write(*,"(A20,F20.16)")'sum_q_initial  ',sum_q_initial
    write(*,"('Relative error:', ES20.4)") (sum_q-sum_q_initial)/sum_q_initial
    write(1110,"(F10.5, ES20.4)") T_temp,(sum_q-sum_q_initial)/sum_q_initial
    
    !!准确解的全局密度通量
    !sum_q_exact = 0.0_prec
    !do i =1,ncells
    !    var_array = 0.0_prec
    !    sum_temp = 0.0_prec
    !    var_array = cellset(i).spvalue_ori_exa(:,:,1)*cellset(i).det_J(:,:)
    !    call Gauss_DInt(var_array,sum_temp)       
    !    sum_q_exact = sum_q_exact + sum_temp
    !end do  
    
    !write(*,"('q_exa',F20.16,ES20.4)") sum_q_exact,abs(sum_q_exact-sum_q_initial)/sum_q_initial
    !write(1110,"(F10.5, F24.16,ES20.4,ES20.4)") T_temp,abs(sum-sum_q_initial)/sum_q_initial,abs(sum_q_exact-sum_q_initial)/sum_q_initial,sum
    !stop
end subroutine solve_errorGCL
subroutine solve_WingFlow_CL
    use real_precision
    use parameter_setting
    use global_var
    use type_module
    use bc_module
    implicit none
    
    integer :: i,j,k,l,m,P1,P2,P3,P4,iL,iR,LC_sideth,RC_sideth,rk
    character(len=100)filename
    character(len=5)char_nsdpx,char_nsdpy,char_nsdx,char_nsdy
    integer,allocatable ::  SC_BP_index(:,:) !SubCell_Bound_Points_Index(4,nsp),取出靠近单元侧边的点
    integer :: Vertex_P_index(4)
    integer :: cells_contain_vertex(10,2),sum_cells_contain_vertex,nextNearSideth,indexNearCell,sideIndex,cellIndex
    integer :: count1,count2,count3,count4,startCell,nextCell,vertex_th,count_temp
    real(prec),dimension(:,:) :: bc_nodes_values(nsp,4)
    integer,dimension(:) :: wallCells_index_set(nbdsides)
    real(prec),dimension(nsp) :: varFunc,varFunc_r
    real(prec) ,external :: Gauss_integral_SPs
    real(prec) :: p_ds,ds,normal_xy(2),r_ave
    
    !!## C_L = L/(0.5*rho*V^2 * S)
    !!## L为升力，V 是远场速度，S是特征面积
    
    !初始化
    CL_xy = 0.0_prec
    !简单寻找出边界附近的求解点
    wallCells_index_set = 0
    count_temp = 1
    do i = 1,nbdsides
        cellIndex = BoundCells_index_set(i)
        if(i>1 .and. cellIndex == BoundCells_index_set(i-1))cycle  !去掉拐角重复的单元
        !write(*,*)cellIndex
        do j = 1,4
            sideIndex = cellset(cellIndex).sides(j)
            if(sideset(sideIndex).bc==Wall)then
                !write(*,*)cellIndex,j
                if(j == 1)then           
                    varFunc = cellset(cellIndex).spvalue_ori(1,:,4)
                    p_ds = Gauss_integral_SPs(varFunc)/cellset(cellIndex).fpMdirect_G(1,1,1)
                    !ds = (xy_coor(cellset(cellIndex).nodes(1),1)-xy_coor(cellset(cellIndex).nodes(2),1))**2+(xy_coor(cellset(cellIndex).nodes(1),2)-xy_coor(cellset(cellIndex).nodes(2),2))**2
                    !ds = sqrt(ds)
                    call ds_normalVec(cellIndex,j,normal_xy(1),normal_xy(2))
                    !varFunc_r = cellset(cellIndex).spvalue_ori(1,:,1)
                    !r_ave = Gauss_integral_SPs(varFunc_r)/2.0_prec
                elseif(j == 2)then
                    varFunc = cellset(cellIndex).spvalue_ori(:,nsp,4)
                    p_ds = Gauss_integral_SPs(varFunc)/cellset(cellIndex).fpMdirect_F(1,nsp+1,4)
                    !ds = (xy_coor(cellset(cellIndex).nodes(2),1)-xy_coor(cellset(cellIndex).nodes(3),1))**2+(xy_coor(cellset(cellIndex).nodes(2),2)-xy_coor(cellset(cellIndex).nodes(3),2))**2
                    !ds = sqrt(ds)
                    call ds_normalVec(cellIndex,j,normal_xy(1),normal_xy(2))                   
                    !varFunc_r = cellset(cellIndex).spvalue_ori(:,nsp,1)
                    !r_ave = Gauss_integral_SPs(varFunc_r)/2.0_prec
                   
                elseif(j == 3)then
                    varFunc = cellset(cellIndex).spvalue_ori(nsp,:,4)
                    p_ds = Gauss_integral_SPs(varFunc)/cellset(cellIndex).fpMdirect_G(1,nsp+1,1)
                    !ds = (xy_coor(cellset(cellIndex).nodes(3),1)-xy_coor(cellset(cellIndex).nodes(4),1))**2+(xy_coor(cellset(cellIndex).nodes(3),2)-xy_coor(cellset(cellIndex).nodes(4),2))**2
                    !ds = sqrt(ds)
                    call ds_normalVec(cellIndex,j,normal_xy(1),normal_xy(2))
                    !varFunc_r = cellset(cellIndex).spvalue_ori(nsp,:,1)
                    !r_ave = Gauss_integral_SPs(varFunc_r)/2.0_prec
                   
                elseif(j == 4)then
                    varFunc = cellset(cellIndex).spvalue_ori(:,1,4)
                    p_ds = Gauss_integral_SPs(varFunc)/cellset(cellIndex).fpMdirect_F(1,1,4)
                    !ds = (xy_coor(cellset(cellIndex).nodes(4),1)-xy_coor(cellset(cellIndex).nodes(1),1))**2+(xy_coor(cellset(cellIndex).nodes(4),2)-xy_coor(cellset(cellIndex).nodes(1),2))**2
                    !ds = sqrt(ds)
                    call ds_normalVec(cellIndex,j,normal_xy(1),normal_xy(2))
                    !varFunc_r = cellset(cellIndex).spvalue_ori(:,1,1)
                    !r_ave = Gauss_integral_SPs(varFunc_r)/2.0_prec                  
                end if               
                CL_xy(:) = CL_xy(:) + p_ds*normal_xy(:)
                !write(*,*) p_ds,p_ds*normal_xy(2)
                exit
            end if
        end do       
    end do
    !write(*,*)'sum', CL_xy
    CL_xy(:) = CL_xy(:)/(0.5_prec*1.4_prec*Ms**2)
    write(1140,*) nt_temp,CL_xy(:)
    !stop
    !write(*,*)CL_xy(:)
    !stop
end subroutine solve_WingFlow_CL
    
subroutine Trouble_dis

    use real_precision
    use parameter_setting
    use global_var
    use type_module
    implicit none
    integer :: i,j,k,l
    !伪一维的情况
    do i = 1, nsdx
        if(cellset(i).Beta == 1)then
            write(1120,*) cellset(i).sp_coor(1,3,1) ,T_temp
        end if
    end do
end subroutine Trouble_dis
    
subroutine Trouble_Proportion   
    use real_precision
    use parameter_setting
    use global_var
    use type_module
    implicit none
    integer :: i,j,k,l,counter_beta,counter_beta_kesi,counter_beta_eta
    real(prec) :: percentage,percentage_kesi,percentage_eta
    !问题单元比例
    counter_beta = 0
    counter_beta_kesi = 0
    counter_beta_eta  = 0
    do i = 1, ncells
        if(cellset(i).Beta==1)then
            counter_beta = counter_beta + 1
        end if
        do k = 1,nsp
            if(cellset(i).Beta_line(k,1) == 1)then
                counter_beta_kesi = counter_beta_kesi + 1
            end if
            if(cellset(i).Beta_line(k,2) == 1)then
                counter_beta_eta = counter_beta_eta + 1
            end if
        end do
    end do
    percentage = real(counter_beta)/real(ncells)
    percentage_kesi = real(counter_beta_kesi)/real(ncells*nsp)
    percentage_eta = real(counter_beta_eta)/real(ncells*nsp)
    write(1121,"(F20.15,3F15.8)") T_temp,percentage,percentage_kesi,percentage_eta
    
end subroutine Trouble_Proportion   
    
subroutine solve_dt

    !-----------------------------------------------------------------------------
    !
    !   方法：求解时间步长
    !   描述：若采取自适应方法，则利用CFL限制求解每一时间步的时间步长，若是固定时间步长，则不进行计算
    !         采取自适应 dt：
    !         dt = CFL * dx / max{|u|+c,|v|+c}
    !         dt = CFL * min{dx/(|u|+c),dy/(|v|+c)}
    !         dt = CFL * min(dx,dy)/max(|u|+c, |v|+c) 采用的
    !   作者：gqShi 
    !   历史：** 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use parameter_setting
    use global_var
    use type_module
    implicit none
    integer :: i
    real(prec) :: c,r,u,v,p,arr_c(nsp,nsp),arr_uc(nsp,nsp),arr_vc(nsp,nsp)
    real(prec) :: dx_local(nsp,nsp),dy_local(nsp,nsp),min_dt
    
    if(dt_switch == dt_adapt)then
        
        min_dt = 1.0_prec
        do i = 1,ncells
            !arr_c = sqrt(gamma*cellset(i).spvalue_ori(:,:,4)/cellset(i).spvalue_ori(:,:,1)) !密度或压力出负就会计算失败
            !arr_uc = abs(cellset(i).spvalue_ori(:,:,2))+arr_c
            !arr_vc = abs(cellset(i).spvalue_ori(:,:,3))+arr_c
            !dx_local = 2.0_prec*cellset(i).MJacobi(:,:,1)
            !dy_local = 2.0_prec*cellset(i).MJacobi(:,:,4)
            !
            !dt = CFL * min(minVal(dx_local),minVal(dy_local))/max(maxVal(arr_uc),maxVal(arr_vc))
            !min_dt = min(min_dt,dt)
            !write(*,*)dt,min(min_dt,dt)
        
            !!! 另一种方法
            !arr_uc = dx_local/arr_uc
            !arr_vc = dy_local/arr_vc
            !
            !dt = CFL * min(minVal(arr_uc),minVal(arr_vc))
            !min_dt = min(min_dt,dt)
            !write(*,*)dt,min(min_dt,dt)
            
            !!!另一种方法
            arr_c = sqrt(gamma*cellset(i).spvalue_ori(:,:,4)/cellset(i).spvalue_ori(:,:,1))
            arr_uc = abs(cellset(i).spvalue_ori(:,:,2))+arr_c
            arr_vc = abs(cellset(i).spvalue_ori(:,:,3))+arr_c
            
            !dt = CFL * min(mindx_global,mindy_global)/max(maxVal(arr_uc),maxVal(arr_vc))
            dt = CFL * min_dis/max(maxVal(arr_uc),maxVal(arr_vc))
            min_dt = min(min_dt,dt)
        end do
        
        dt = min_dt/(4.0_prec*(nsp-1)+1)  !1/d * 1/(2N+1)
    elseif(dt_switch == dt_fixed)then
        
        ! nothing 文件输入dt

    end if
    
end subroutine solve_dt 
    
subroutine ds_normalVec(indexCell,sideth,normalX,normalY)
    !计算(x1,y1)(x2,y2)线段的法向单位向量
    use real_precision
    use global_var
    use type_module
    implicit none
    real(prec) :: x1,y1,x2,y2,xa,ya,x,y,normalX,normalY
    integer :: indexCell,sideth,vertexth1,vertexth2
    integer :: k
    !计算外法向向量
    vertexth1 = sideth
    vertexth2 = sideth+1
    if(vertexth1 == 4) vertexth2 = 1
    x1 = xy_coor(cellset(indexCell).nodes(vertexth1),1)
    y1 = xy_coor(cellset(indexCell).nodes(vertexth1),2)
    x2 = xy_coor(cellset(indexCell).nodes(vertexth2),1)
    y2 = xy_coor(cellset(indexCell).nodes(vertexth2),2)
    do k = 1,4
        if(k .NE. vertexth1 .and. k .NE.vertexth2)then
            xa = xy_coor(cellset(indexCell).nodes(k),1)
            ya = xy_coor(cellset(indexCell).nodes(k),2)
            exit
        end if       
    end do   
    call sym_Point(x1,y1,x2,y2,xa,ya,x,y) 

    normalX = (x-xa)/sqrt((x-xa)**2 + (y-ya)**2)
    normalY = (y-ya)/sqrt((x-xa)**2 + (y-ya)**2)

  
end subroutine ds_normalVec