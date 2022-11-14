subroutine time_advancing

    !-----------------------------------------------------------------------------
    !
    !   ������ʱ���ƽ�ģ��
    !   �������漰�����ݽ������߽紦��ʱ����ɢ��������׽��
    !   ���ߣ�gqShi 
    !   ��ʷ��** 
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
    
    call allocateSpace              !----�����ڴ�  
    call createFile                 !----�����ļ�  
    call update_Boundary_fluent     !----���±߽�
     
    call ori_to_con                 !----ԭʼ����->�غ����
    call ori_to_flu                 !----ԭʼ����->����ͨ��
    
    nt_temp = 1
    T_temp = 0.0_prec
    DetectionTime = 0.0_prec
    
    do while(T_temp < T)
        
        ! ����ʱ�䲽��
        call solve_dt
        T_temp = T_temp + dt
        
        ! �ж���ֹ����
        if(T_temp > T )then  
            dt =  T - (T_temp - dt)
            T_temp = T
            write(*,*)'dt',dt
        end if
        
        ! �в������Ҫ _old 
        do i = 1,ncells
            cellset(i).spvalue_ori_old  =  cellset(i).spvalue_ori 
        end do
        
        ! ��ʱ���ƣ�20��
        if(nt_temp == 1 .OR.nt_temp == 21)then     
            call com_time(nt_temp)
            nstep = T/dt
            write(*,*) 'nstep,dt',nstep,dt
        end if

        
        ! Runge-Kutta
        do RKstage = 1,3        

            call mark_TroCell               !----������ⵥԪ
            !call Trouble_dis               !----αһά������ⵥԪ�ֲ�/��ά������ⵥԪռ��     
            call update_Boundary_Detect     !----���±߽����ⵥԪ����� 
            call phy_to_com                 !----�任������ռ�   
            call Face_Flux_upw              !----����߽�ͨ��   
            
            call RK_advancing(RKstage)      !----RK �ƽ�         
            
            !call one_dimensionalized       !----һά��.�Ų��ά��һά֮��ı仯           
            call con_to_ori                 !----�غ����->ԭʼ����
            call update_Boundary_fluent     !----���±߽�     
            !call correct_density_pressure  !----����ܶȣ�ѹ���Ƿ���ָ�ֵ����������                                    
            call ori_to_flu                 !----ԭʼ����->����ͨ�� 

        end do
        
        !��ӡ�м����ݣ��в�غ�����ֵ���
        call printResultPerStep
        
        !�ƽ�����
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
    !   ��������������������Ҫ���ļ�
    !   ����������ɢ�غ�������ļ������ⵥԪ�ֲ��ļ����в��ļ���������������ϵ�������ļ���
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use preParameter
    implicit none

    ! ������ɢ�غ��������غ���
    if(compuConsLaw == consLaw_open)then
        open(1110,file = 'error_of_global_con_law.plt')    
        sum_q_initial = 0.0_prec
        call exact_solution(0.0_prec)
        call solve_errorGCL(0) 
    endif

    ! ���ⵥԪ��¼,αһά��¼�ֲ���ʱ��仯����ά��¼���ⵥԪռ����ʱ��仯
    if(trouCell_record == trouCell_yes)then
        open(1120,file = 'troubleCellRecord.plt')
        open(1121,file = 'troubleCellProportion.plt')
    endif
    
    ! ����в�
    if(compuResidua == residua_yes)then
        open(1130,file = 'residualsRecord.plt')
        write(1130,*)'Title="residuals"'
        write(1130,*)'Variables= nstep,resRelative,res,resAve'
    endif
    
    ! ������������������ϵ������
    if(case_comp == WingFlow_case)then
        open(1140,file = 'CL_WingFlow.plt')
        CL_xy(:) = 0.0_prec
    end if

end subroutine createFile    
    
subroutine allocateSpace
 
    !-----------------------------------------------------------------------------
    !
    !   ��������������Ҫ�õ�����������ڴ�
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
    !   ��������ӡ�����ļ�
    !   ������ʱ���ƽ������У�ÿ��n����ӡһ����������
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use preParameter
    implicit none
    
    !!��ӡ���ⵥԪ
    if(trouCell_record == trouCell_yes)then
        
        ! αһά�������y������ݵ�һ������ֵ
        !call Trouble_dis
        
        ! ��ά������ⵥԪռ��
        call Trouble_Proportion
    end if
            
    ! �в�
    if(compuResidua == residua_yes)then
        if(nt_temp == 1 .or. mod(nt_temp,20)==0)then    !----ÿ20�����һ����ɢ�غ����
            call solve_residuals(nt_temp)
        end if
    endif
    
    ! ��ɢ�غ����
    if(compuConsLaw == consLaw_open .and. mod(nt_temp,20)==0 )then
        call solve_errorGCL(nt_temp)                    !----ÿ20�����һ����ɢ�غ����
    endif
        
    ! ��������������ϵ������
    if(case_comp == WingFlow_case .and. mod(nt_temp,20)==0)then
        call solve_WingFlow_CL                          !----ÿ20�����һ����ɢ�غ����
    end if
        
    ! ������
    if(mod(nt_temp,print_step)==0)then
                
        write(*,"('step: 'I8',   time'F16.5',   dt'E16.5)") nt_temp,T_temp,dt
        call now_time                                   !----����nt_temp����ʱ��
        call print_num_data(T_temp)                     !----ÿprint_step�� ���һ����ֵ��
        !call exact_solution(T_temp)                    !----���׼ȷ��
        !call print_exa_data(T_temp)                    !----��ӡ׼ȷ��
        !call renewal_program                           !----�����ļ���δ�ӡ�û���ü����������           
    end if
    
end subroutine printResultPerStep
    
subroutine solve_residuals(nt_t)
 
    !-----------------------------------------------------------------------------
    !
    !   �������в����
    !   ����������֮����ܶ������Ե�һ������������
    !   ���ߣ�gqShi 
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
    !   ������Runge-Kuttaʱ���ƽ�
    !   ����������ʵ������׵�RK�������ڴ�������TVD RK
    !   ���ߣ�gqShi 
    !   ��ʷ��** //2021.12.08
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
            ! �������Ƿ�Ϊ��
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
    !   ��������άһά��
    !   ���������ṹ����αһά����ת����ÿ��һά��Ԫ����ȣ��Ų��ά���ص�Ӱ�죬���Լӵģ��ò������ù�
    !   ���ߣ�gqShi 
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
        !x ����
        do i = 1,nsdx
            do j = 1,nsdy
                do k = 1,nsp
                    cellset(i+(j-1)*nsdx).spvalue_con(k,:,:) = cellset(i).spvalue_con(1,:,:)
                end do        
            end do
        end do
    case(direct_y)
        !y ����
        do i = 1,nsdx
            do j = 1,nsdy
                do k = 1,nsp
                    cellset(i+(j-1)*nsdx).spvalue_con(:,k,:) = cellset(1+(j-1)*nsdx).spvalue_con(:,1,:)
                end do        
            end do
        end do
    end select

end subroutine one_dimensionalized