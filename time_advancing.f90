subroutine time_advancing

    !-----------------------------------------------------------------------------
    !
    !   方法：时间推进模块
    !   描述：涉及到数据交换，边界处理，时间离散，激波捕捉等
    !   作者：gqShi 
    !   历史：** 
    !
    !-----------------------------------------------------------------------------
    
    use parameter_setting
    use global_var
    use type_module
    use time_test
    use preParameter
    use Time_module
    implicit none
    
    integer nt,i,j,k,l,RKstage
    integer :: int_time1,int_time2,int_time3,int_time4,int_time5,int_time6,int_time7

    write(*,*) 'Time advance.'
    
    call allocateSpace              !----分配内存  
    call createFile                 !----创建文件  
    call update_Boundary_fluent     !----更新边界
     
    call ori_to_con                 !----原始变量->守恒变量
    call ori_to_flu                 !----原始变量->对流通量
    
    nt_temp = 1
    T_temp = 0.0_prec
    DetectionTime = 0.0_prec
    
    do while(T_temp < T)
        
        ! 更新时间步长
        call solve_dt
        T_temp = T_temp + dt
        
        ! 判断终止计算
        if(T_temp > T )then  
            dt =  T - (T_temp - dt)
            T_temp = T
            write(*,*)'dt',dt
        end if
        
        ! 残差计算需要 _old 
        do i = 1,ncells
            cellset(i).spvalue_ori_old  =  cellset(i).spvalue_ori 
        end do
        
        ! 耗时估计，20步
        if(nt_temp == 1 .OR.nt_temp == 21)then     
            call com_time(nt_temp)
            nstep = T/dt
            write(*,*) 'nstep,dt',nstep,dt
        end if

        
        ! Runge-Kutta
        do RKstage = 1,3        

            call mark_TroCell               !----侦测问题单元
            !call Trouble_dis               !----伪一维情况问题单元分布/二维情况问题单元占比     
            call update_Boundary_Detect     !----更新边界虚拟单元的侦测 
            call phy_to_com                 !----变换到计算空间   
            call Face_Flux_upw              !----计算边界通量   
            
            call RK_advancing(RKstage)      !----RK 推进         
            
            !call one_dimensionalized       !----一维化.排查二维与一维之间的变化           
            call con_to_ori                 !----守恒变量->原始变量
            call update_Boundary_fluent     !----更新边界     
            !call correct_density_pressure  !----检查密度，压力是否出现负值。矫正程序                                    
            call ori_to_flu                 !----原始变量->对流通量 

        end do
        
        !打印中间数据，残差，守恒误差，数值解等
        call printResultPerStep
        
        !推进步数
        nt_temp = nt_temp + 1
    end do
    
    write(*,*)'Final step:',nt_temp,T_temp,dt
    
    close(1110)
    close(1120)
    close(1121)
    close(1130)
    close(1140)
    
end subroutine time_advancing
    
subroutine createFile    
 
    !-----------------------------------------------------------------------------
    !
    !   方法：创建求解过程中需要的文件
    !   描述：如离散守恒律误差文件、问题单元分布文件、残差文件、翼型绕流升力系数数据文件等
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use preParameter
    implicit none

    ! 计算离散守恒误差，测试守恒性
    if(compuConsLaw == consLaw_open)then
        open(1110,file = 'error_of_global_con_law.plt')    
        sum_q_initial = 0.0_prec
        call exact_solution(0.0_prec)
        call solve_errorGCL(0) 
    endif

    ! 问题单元记录,伪一维记录分布随时间变化，二维记录问题单元占比随时间变化
    if(trouCell_record == trouCell_yes)then
        open(1120,file = 'troubleCellRecord.plt')
        open(1121,file = 'troubleCellProportion.plt')
    endif
    
    ! 计算残差
    if(compuResidua == residua_yes)then
        open(1130,file = 'residualsRecord.plt')
        write(1130,*)'Title="residuals"'
        write(1130,*)'Variables= nstep,resRelative,res,resAve'
    endif
    
    ! 计算翼型绕流，升力系数曲线
    if(case_comp == WingFlow_case)then
        open(1140,file = 'CL_WingFlow.plt')
        CL_xy(:) = 0.0_prec
    end if

end subroutine createFile    
    
subroutine allocateSpace
 
    !-----------------------------------------------------------------------------
    !
    !   方法：给接下来要用到的数组分配内存
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use preParameter
    implicit none
    
    integer i

   
    do i = 1,ncells+nbdsides
        allocate(cellset(i).spvalue_con(nsp,nsp,4),cellset(i).spvalue_fluF(nsp,nsp,4),cellset(i).spvalue_fluG(nsp,nsp,4))
        allocate(cellset(i).spvalue_con_loc(nsp,nsp,4),cellset(i).spvalue_fluF_loc(nsp,nsp,4),cellset(i).spvalue_fluG_loc(nsp,nsp,4))
        allocate(cellset(i).spvalue_con_tem(nsp,nsp,4))
        allocate(cellset(i).fluxF_innerfp(nsp,nsp+1,4),cellset(i).fluxG_innerfp(nsp,nsp+1,4))
        allocate(cellset(i).spvalue_ori_exa(nsp,nsp,4))
        allocate(cellset(i).spvalue_ori_old(nsp,nsp,4))
        allocate(cellset(i).Beta_line(nsp,2))
        allocate(cellset(i).Smooth_line_x(nsp,nsp+1),cellset(i).Smooth_line_y(nsp,nsp+1))
        cellset(i).Beta_line = 0
    end do
    
    do i = 1,nsides
        allocate(sideset(i).fpvalue_upw(nsp,4))
        sideset(i).fpvalue_upw(:,:) = 0.0_prec 
    end do
    
end subroutine allocateSpace
    
subroutine printResultPerStep
 
    !-----------------------------------------------------------------------------
    !
    !   方法：打印过程文件
    !   描述：时间推进过程中，每隔n步打印一次流场数据
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use preParameter
    implicit none
    
    !!打印问题单元
    if(trouCell_record == trouCell_yes)then
        
        ! 伪一维情况。如y方向根据第一层结果赋值
        !call Trouble_dis
        
        ! 二维输出问题单元占比
        call Trouble_Proportion
    end if
            
    ! 残差
    if(compuResidua == residua_yes)then
        if(nt_temp == 1 .or. mod(nt_temp,20)==0)then    !----每20步输出一次离散守恒误差
            call solve_residuals(nt_temp)
        end if
    endif
    
    ! 离散守恒误差
    if(compuConsLaw == consLaw_open .and. mod(nt_temp,20)==0 )then
        call solve_errorGCL(nt_temp)                    !----每20步输出一次离散守恒误差
    endif
        
    ! 翼型绕流，升力系数曲线
    if(case_comp == WingFlow_case .and. mod(nt_temp,20)==0)then
        call solve_WingFlow_CL                          !----每20步输出一次离散守恒误差
    end if
        
    ! 结果输出
    if(mod(nt_temp,print_step)==0)then
                
        write(*,"('step: 'I8',   time'F16.5',   dt'E16.5)") nt_temp,T_temp,dt
        call now_time                                   !----运行nt_temp步的时间
        call print_num_data(T_temp)                     !----每print_step步 输出一次数值解
        !call exact_solution(T_temp)                    !----求解准确解
        !call print_exa_data(T_temp)                    !----打印准确解
        !call renewal_program                           !----续算文件，未加。没来得及加续算程序           
    end if
    
end subroutine printResultPerStep
    
subroutine solve_residuals(nt_t)
 
    !-----------------------------------------------------------------------------
    !
    !   方法：残差计算
    !   描述：两步之间的密度误差相对第一步误差的相对误差
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    implicit none
    
    integer :: i,j,k,l,nt_t
    real(prec) :: maxResiduals_temp,maxResiduals_Relative,aveResidual
    
    maxResiduals_temp = 0.0_prec
    maxResiduals = 0.0_prec
    maxResiduals_Relative = 0.0_prec
    aveResidual = 0.0_prec
    do i = 1,ncells
        maxResiduals_temp = maxVal(abs(cellset(i).spvalue_ori_old(:,:,1) - cellset(i).spvalue_ori(:,:,1))) 
        maxResiduals = max(maxResiduals,maxResiduals_temp)
        aveResidual = aveResidual + sum(cellset(i).spvalue_ori_old(:,:,1) - cellset(i).spvalue_ori(:,:,1))
    end do
    if(nt_t == 1)then
        FirstResi = maxResiduals
        FirstAveResi = aveResidual/ncells/nsp**2
    end if
    
    if(mod(nt_t,20)==0 .or. nt_t == 1)then
        maxResiduals_Relative = maxResiduals/(FirstResi+1.0e-16)
        aveResidual = aveResidual/ncells/nsp**2
        write(*,"('Residuals: '3F24.16)") maxResiduals_Relative,maxResiduals,aveResidual/FirstAveResi
        write(1130,*) nt_t,maxResiduals_Relative,maxResiduals,aveResidual/FirstAveResi!,T_temp
    end if
    
end subroutine solve_residuals

subroutine RK_advancing(stage)
     
    !-----------------------------------------------------------------------------
    !
    !   方法：Runge-Kutta时间推进
    !   描述：可以实现任意阶的RK方法，在此是三阶TVD RK
    !   作者：gqShi 
    !   历史：** //2021.12.08
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use global_var
    use type_module
    implicit none
    integer i,stage
    integer,parameter :: stage1 = 1,stage2 = 2,stage3 = 3
    real(prec),dimension(:,:,:) :: Lsub(nsp,nsp,4)
    Lsub = 0.0_prec
    
    select case(stage)
    case(stage1)
        do i = 1,ncells
            ! RK1
            call RHS_CPR_NNW(i,Lsub)
            cellset(i).spvalue_con_tem(:,:,:) = cellset(i).spvalue_con(:,:,:)
            cellset(i).spvalue_con(:,:,:) = cellset(i).spvalue_con_tem(:,:,:) + dt*Lsub(:,:,:)
            ! 计算结果是否为正
            if(isnan(cellset(i).spvalue_con(1,1,1)))then
                write(*,"('error! NAN occur in',I6,'-th step,',I6,'-th cell')")nt_temp,i
                call print_num_data(T_temp)
                stop
            end if
        end do
    case(stage2)
        do i = 1,ncells
            ! RK2
            call RHS_CPR_NNW(i,Lsub)
            cellset(i).spvalue_con(:,:,:) = 3.0_prec/4.0_prec*cellset(i).spvalue_con_tem(:,:,:) + 1.0_prec/4.0_prec*(cellset(i).spvalue_con(:,:,:)+dt*Lsub(:,:,:))
        end do        
    case(stage3)
        do i = 1,ncells
            ! RK3
            call RHS_CPR_NNW(i,Lsub)
            cellset(i).spvalue_con(:,:,:) = 1.0_prec/3.0_prec*cellset(i).spvalue_con_tem(:,:,:) + 2.0_prec/3.0_prec*(cellset(i).spvalue_con(:,:,:)+dt*Lsub(:,:,:))
        end do
    case default
        write(*,*)'Error! Occur in subroutine RK_advancing.'
        stop
    end select
    
end subroutine RK_advancing
    
subroutine one_dimensionalized
 
    !-----------------------------------------------------------------------------
    !
    !   方法：二维一维化
    !   描述：将结构网格，伪一维算例转化成每个一维单元都相等，排查二维因素的影响，调试加的，用不到不用管
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer :: i,j,k,l
      
    select case(dire_shock)
    case(direct_x)
        !x 方向
        do i = 1,nsdx
            do j = 1,nsdy
                do k = 1,nsp
                    cellset(i+(j-1)*nsdx).spvalue_con(k,:,:) = cellset(i).spvalue_con(1,:,:)
                end do        
            end do
        end do
    case(direct_y)
        !y 方向
        do i = 1,nsdx
            do j = 1,nsdy
                do k = 1,nsp
                    cellset(i+(j-1)*nsdx).spvalue_con(:,k,:) = cellset(1+(j-1)*nsdx).spvalue_con(:,1,:)
                end do        
            end do
        end do
    end select

end subroutine one_dimensionalized