subroutine mark_TroCell
     
    !-----------------------------------------------------------------------------
    !
    !   ��������Ⲣ������ⵥԪ
    !   ������
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use Time_module
    implicit none
    
    integer :: i,j,k,l

    !call cpu_time(DetectionTimeStart)           !----��⻨��ʱ��. START
    
    !-----------------------------------------------------------------------------
    !
    !   ���������ݼ�������ѡ���Ƿ�������⺯��
    !
    !-----------------------------------------------------------------------------
    
    if(scheme_kind == scheme_cpr)then
        
        do i = 1,ncells+nbdsides
            cellset(i).Beta = 0
            cellset(i).Beta_line = 0
            cellset(i).Smooth_line_x = 0
            cellset(i).Smooth_line_y = 0
        end do
        
    elseif(scheme_kind == scheme_two)then
        
        do i = 1,ncells+nbdsides
            cellset(i).Beta = 1
            cellset(i).Beta_line = 1
            cellset(i).Smooth_line_x = 1
            cellset(i).Smooth_line_y = 1
        end do
        
    elseif(scheme_kind == scheme_hybrid)then
        
        do i = 1,ncells+nbdsides
            cellset(i).Beta_old = cellset(i).Beta
            cellset(i).Beta = 0
            cellset(i).Beta_line = 0
        end do
        
        call cpu_time(DetectionTimeStart)           !----��⻨��ʱ��. START
        
        ! ѡ�������
        select case(detection_type)
        case(detection_TVB)      
            call indicator_TVB              
        case(detection_ATV)          
            call indicator_ATV
        case(detection_MV)       
            call indicator_MV
        case(detection_MDH)
            call indicator_MDH
        case(detection_MDHm)
            call indicator_MDH3
        case(detection_KXRCF)
            call indicator_KXRCF
        case(detection_JST)
            call indicator_JST
        case default
            ! call indicator_test
        end select
    
        call cpu_time(DetectionTimeEnd)             !----��⻨��ʱ��. END
        DetectionTime = DetectionTime + (DetectionTimeEnd-DetectionTimeStart)
        
        ! �Ƿ���ӹ��ɵ�Ԫ
        if(buffer_cell_switch == buffer_cell_yes)then
            call add_buffer_cell
        end if
        
        ! Ȩֵ��������
        call post_weight_func
        
    end if
    
    
    
end subroutine mark_TroCell

subroutine add_buffer_cell
 
    !-----------------------------------------------------------------------------
    !
    !   ��������ӹ��ɵ�Ԫ
    !   �����������ⵥԪ�Ĺ������ڵ�ԪҲ����Ϊ���ⵥԪ
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    
    integer,dimension(:),allocatable :: Cell_Beta
    integer,dimension(:,:,:),allocatable :: Cell_Beta_line
    integer :: i,j,k,l,m,indexNearCell
    integer :: sideth(4),near_kk(4)
    
    allocate(Cell_Beta(ncells),Cell_Beta_line(ncells,nsp,2))
    do i = 1,ncells
        Cell_Beta(i) = cellset(i).Beta
        Cell_Beta_line(i,:,:) = cellset(i).Beta_line(:,:)
    end do
    
    if(detect_type == ByCell)then
        do i = 1,ncells
            if(Cell_Beta(i)==0)cycle    !�����������ⵥԪ��ֱ������
            do j = 1,4
                cellset(cellset(i).nearcells(j)).Beta = 1   
                cellset(cellset(i).nearcells(j)).Beta_line = 1  
            end do
        end do
    elseif(detect_type == ByDim)then
        do i = 1,ncells         
            if(Cell_Beta(i)==0)cycle    !�����������ⵥԪ��ֱ������
            
            !-----��������i�����ⵥԪ�����         
            do j = 1,4  !4�����ڵ�Ԫ
                indexNearCell  = cellset(i).nearcells(j)!���ڵ�Ԫ������
                do k = 1,4
                    if(i==cellset(indexNearCell).nearcells(k))then
                        sideth(j) = k   !��¼����Ԫ�Ĳ�������ڵ�Ԫ �ĵڼ����
                    end if
                end do        
            end do 
            
            do k =1,nsp
                do j = 1,4
                    if(sideth(j)*j==2.or.sideth(j)*j==12.or.sideth(j)==j)then
                        near_kk(j) = nsp+1-k    !�ж��ڵ�Ԫ�ǵڼ���
                    else 
                        near_kk(j) = k
                    end if
                end do
                
                !x����
                if(Cell_Beta_line(i,k,1) == 1)then                                      
                    if(sideth(4)==2.or.sideth(4)==4)then
                        cellset(cellset(i).nearcells(4)).Beta_line(near_kk(4),1) = 1                       
                    elseif(sideth(4)==1.or.sideth(4)==3)then
                        cellset(cellset(i).nearcells(4)).Beta_line(near_kk(4),2) = 1
                    end if
                    if(sideth(2)==2.or.sideth(2)==4)then
                        cellset(cellset(i).nearcells(2)).Beta_line(near_kk(2),1) = 1
                    elseif(sideth(2)==1.or.sideth(2)==3)then
                        cellset(cellset(i).nearcells(2)).Beta_line(near_kk(2),2) = 1
                    end if
                    cellset(cellset(i).nearcells(2)).Beta = 1
                    cellset(cellset(i).nearcells(4)).Beta = 1
                end if
                
                !y����
                if(Cell_Beta_line(i,k,2) == 1)then
                    if(sideth(1)==2.or.sideth(1)==4)then
                        cellset(cellset(i).nearcells(1)).Beta_line(near_kk(1),1) = 1
                    elseif(sideth(1)==1.or.sideth(1)==3)then
                        cellset(cellset(i).nearcells(1)).Beta_line(near_kk(1),2) = 1
                    end if
                    if(sideth(3)==2.or.sideth(3)==4)then
                        cellset(cellset(i).nearcells(3)).Beta_line(near_kk(3),1) = 1
                    elseif(sideth(3)==1.or.sideth(3)==3)then
                        cellset(cellset(i).nearcells(3)).Beta_line(near_kk(3),2) = 1
                    end if     
                    cellset(cellset(i).nearcells(1)).Beta = 1
                    cellset(cellset(i).nearcells(3)).Beta = 1
                end if                
            end do
        end do   
    end if
   
    deallocate(Cell_Beta,Cell_Beta_line)
    
end subroutine add_buffer_cell

subroutine post_weight_func
 
    !-----------------------------------------------------------------------------
    !
    !   ������
    !   ��������subroutine get_oriL_weight���á�
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    
    integer ::index,nearCellIndex,LC_sideth,RC_sideth,l,m,k,RC_l
    real(prec),dimension(:,:) :: u(3,4),conL(4),ll(4,4),rr(4,4),chara_con(3,4),ui_all(8,4)
    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4),varConL(nsp,nsp,4),varConR(nsp,nsp,4)
    real(prec),external :: LaI_nPs
    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),oriL(4),chara_oriL(4),chara_conL(4),uA(4)
    real(prec) :: realDis1,realDis2,detJ_fp,dP1P2
    integer :: beta_dim,beta_dim_near
    real(prec),dimension(:) :: oriL_CPR(4),oriL_CNNW2(4)
    real(prec) :: smoother
    integer :: i,j
    
    do i = 1,ncells
        do l = 1,nsp  
            !x 
            if(cellset(i).Beta_line(l,1) == 1)then
                cycle
            else
                !4
                LC_sideth = 4
                nearCellIndex = cellset(i).nearcells(LC_sideth)
                do j = 1,4                 
                    if(cellset(nearCellIndex).nearcells(j) == i .and. nearCellIndex < ncells)then
                        RC_sideth = j
                        exit
                    end if
                end do
                call get_RC_l(LC_sideth,RC_sideth,l,RC_l)
                if(RC_sideth == 2 .or.RC_sideth == 4)then  
                    beta_dim_near = cellset(nearCellIndex).Beta_line(RC_l,1)
                elseif(RC_sideth == 1 .or.RC_sideth == 3)then
                    beta_dim_near  = cellset(nearCellIndex).Beta_line(RC_l,2)              
                end if
                if(beta_dim_near == 1)then
                    cycle
                endif
                !2
                LC_sideth = 2
                nearCellIndex = cellset(i).nearcells(LC_sideth)
                do j = 1,4                 
                    if(cellset(nearCellIndex).nearcells(j) == i .and. nearCellIndex < ncells)then
                        RC_sideth = j
                        exit
                    end if
                end do
                call get_RC_l(LC_sideth,RC_sideth,l,RC_l)
                if(RC_sideth == 2 .or.RC_sideth == 4)then  
                    beta_dim_near = cellset(nearCellIndex).Beta_line(RC_l,1)
                elseif(RC_sideth == 1 .or.RC_sideth == 3)then
                    beta_dim_near  = cellset(nearCellIndex).Beta_line(RC_l,2)              
                end if
                if(beta_dim_near == 1)then
                    cycle
                else
                    cellset(i).Smooth_line_x(l,:) = 0
                endif
            endif                        
        end do
        
        do l = 1,nsp  
          !y 
            if(cellset(i).Beta_line(l,2) == 1)then
                cycle
            else
                !1
                LC_sideth = 1
                nearCellIndex = cellset(i).nearcells(LC_sideth)
                do j = 1,4                 
                    if(cellset(nearCellIndex).nearcells(j) == i .and. nearCellIndex < ncells)then
                        RC_sideth = j
                        exit
                    end if
                end do
                call get_RC_l(LC_sideth,RC_sideth,l,RC_l)
                if(RC_sideth == 2 .or.RC_sideth == 4)then  
                    beta_dim_near = cellset(nearCellIndex).Beta_line(RC_l,1)
                elseif(RC_sideth == 1 .or.RC_sideth == 3)then
                    beta_dim_near  = cellset(nearCellIndex).Beta_line(RC_l,2)              
                end if
                if(beta_dim_near == 1)then
                    cycle
                endif
                !3
                LC_sideth = 3
                nearCellIndex = cellset(i).nearcells(LC_sideth)
                do j = 1,4                 
                    if(cellset(nearCellIndex).nearcells(j) == i .and. nearCellIndex < ncells)then
                        RC_sideth = j
                        exit
                    end if
                end do
                call get_RC_l(LC_sideth,RC_sideth,l,RC_l)
                if(RC_sideth == 2 .or.RC_sideth == 4)then  
                    beta_dim_near = cellset(nearCellIndex).Beta_line(RC_l,1)
                elseif(RC_sideth == 1 .or.RC_sideth == 3)then
                    beta_dim_near  = cellset(nearCellIndex).Beta_line(RC_l,2)              
                end if
                if(beta_dim_near == 1)then
                    cycle
                else
                    cellset(i).Smooth_line_y(l,:) = 0
                endif
            endif                         
        end do
    end do
    
end subroutine post_weight_func