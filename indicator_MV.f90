subroutine indicator_MV
    !Refs:[1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29.
    !   ʵ���ϣ���ATV����С�ĸĶ���ԭ�����ǲ�ȡAverage Total Variation,�������������ء��ڴ˽�ƽ���ܱ���Ϊ����
    !   ������ˣ�Ҳ�������⣬�Ǿ��ǹ⻬����������ⵥԪ������ǿ�������ڵ�����£����ܶ��������ⲻ������
    !   ������֤������1����ȡȫ��������ֵ*theta��Ϊ�ж����ݣ�2����ȡ��ȫ��������ֵ-��Сֵ��*theta��Ϊ�ж����ݡ�
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,near_k,ikesi1,ikesi2,ieta1,ieta2,ikesi,ieta
    integer :: sideth(4),near_kk(4),sidethL,sidethR,indexNearCell
    integer :: index,indexCellL,indexCellR,Beta
    real(prec),dimension(4) :: Max_TV,Max_v,Min_v
    real(prec),dimension(nsp,nsp,4) :: varOriCell,varOriCellL,varOriCellR
    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
    real(prec) :: temp,theta,Vtemp_x(4),Vtemp_y(4)
    
    theta = 0.5_prec
    ! �ȼ���TV��the total variation of the solution��.��������֮�䶼��Ҫ���㡣
    Max_TV = 0.0_prec
    Max_v = 0.0_prec
    Min_v = 0.0_prec
    !���ȼ���һ����Ԫ�ڲ������TV
    do i = 1,ncells               
        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)!!����ȡ����Ԫ -- ������ԭʼ����             
        do m = 1,1  !������
            Max_v = MaxVal(varOriCell(:,:,m))
            Min_v = MinVal(varOriCell(:,:,m))
            do k = 1,nsp
                do l = 1,nsp-1
                    temp = abs(varOriCell(k,l+1,m) - varOriCell(k,l,m)) !x direction
                    Max_TV(m) = max(Max_TV(m),temp)
                    temp = abs(varOriCell(l+1,k,m) - varOriCell(l,k,m)) !y
                    Max_TV(m) = max(Max_TV(m),temp)
                end do
            end do      
        end do
    end do        
    
    !Ȼ����ݱ߼���ߵ����ڵ�Ԫ�ڽ��ӱߵ��������
    do i = 1,nsides        
        indexCellL = sideset(i).nearcells(1)
        indexCellR = sideset(i).nearcells(2)
        if(indexCellL == 0 .OR. indexCellR == 0)cycle
        varOriCellL(:,:,:) = cellset(indexCellL).spvalue_ori(:,:,:)
        varOriCellR(:,:,:) = cellset(indexCellR).spvalue_ori(:,:,:)
        do j = 1,4
            if(i == cellset(indexCellL).sides(j)) sidethL = j            !���������ڵ�Ԫ�ĵڼ����
            if(i == cellset(indexCellR).sides(j)) sidethR = j
        end do
                
        do m = 1,1
            do k = 1,nsp
                if((sidethL*sidethR==2).or.(sidethL*sidethR==12).or.(sidethL==sidethR))then!������������������������෴
                    near_k = nsp+1-k
                else
                    near_k = k
                end if
                if(sidethL==1)then
                    ikesi1 = 1;     ieta1 = k;
                elseif(sidethL==2)then
                    ikesi1 = k;     ieta1 = nsp
                elseif(sidethL==3)then
                    ikesi1 = nsp;   ieta1 = k;
                elseif(sidethL==4)then
                    ikesi1 = k;     ieta1 = 1
                end if
                if(sidethR==1)then
                    ikesi2 = 1;     ieta2 = near_k
                elseif(sidethR==2)then
                    ikesi2 = near_k;ieta2 = nsp
                elseif(sidethR==3)then
                    ikesi2 = nsp;   ieta2 = near_k
                elseif(sidethR==4)then
                    ikesi2 = near_k;ieta2 = 1
                end if
                temp = abs(varOriCellL(ikesi1,ieta1,m)-varOriCellR(ikesi2,ieta2,m))               
                Max_TV(m) = max(Max_TV(m),temp)
            end do
        end do
    end do

    !-------------------------------------------------------
    !!����Ԫ -- ������ԭʼ����
    do i = 1,ncells
        Vtemp_x = 0.0_prec
        Vtemp_y = 0.0_prec
        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)    
        !�ڵ�Ԫ
        do j = 1,4  !4�����ڵ�Ԫ
            indexNearCell  = cellset(i).nearcells(j)!���ڵ�Ԫ������
            do k = 1,4
                if(i==cellset(indexNearCell).nearcells(k))then
                    sideth(j) = k   !��¼����Ԫ�Ĳ�������ڵ�Ԫ �ĵڼ����
                end if
            end do 
            varOriCellNearCell(j,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
        end do 
        !----x direction --------------------------------------------------------------------
        beta = 0
        do m = 1,1
            !write(*,*)'12'
            do k =1,nsp
                if(sideth(4)==1)then    !x���򣬿������4�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
                    ikesi = 1;          ieta = k;
                elseif(sideth(4)==2)then
                    ikesi = k;    ieta = nsp
                elseif(sideth(4)==3)then
                    ikesi = nsp;   ieta = nsp+1-k;
                elseif(sideth(4)==4)then
                    ikesi = nsp+1-k;     ieta = 1
                end if  
   
                temp =  abs(varOriCell(k,1,m)-varOriCellNearCell(4,ikesi,ieta,m))
                Vtemp_x(m) = max(Vtemp_x(m),temp)                
                
                do l = 2,nsp
                    temp =  abs(varOriCell(k,l,m)-varOriCell(k,l-1,m))
                    Vtemp_x(m) = max(Vtemp_x(m),temp)
                    if(temp>theta*Max_TV(m))then    !MV : theta*Max_TV(m);  
                        cellset(i).Beta_line(k,1) = 1
                        cellset(i).Beta = 1
                    end if
                end do 
                !x���򣬿������2�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
                if(sideth(2)==1)then    
                    ikesi = 1;          ieta = nsp+1-k;
                elseif(sideth(2)==2)then
                    ikesi = nsp+1-k;    ieta = nsp
                elseif(sideth(2)==3)then
                    ikesi = nsp;   ieta = k;
                elseif(sideth(2)==4)then
                    ikesi = k;     ieta = 1
                end if  
                temp =  abs(varOriCell(k,nsp,m)-varOriCellNearCell(2,ikesi,ieta,m))
                Vtemp_x(m) = max(Vtemp_x(m),temp)

            !--------------------------------------------------------------------------------------------------       
    
            !----y direction ----------------------------------------------------------------------------------
                if(sideth(1)==1)then    !y���򣬿������1�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
                    ikesi = 1;          ieta = nsp+1-k;
                elseif(sideth(1)==2)then
                    ikesi = nsp+1-k;    ieta = nsp
                elseif(sideth(1)==3)then
                    ikesi = nsp;   ieta = k;
                elseif(sideth(1)==4)then
                    ikesi = k;     ieta = 1
                end if  
                temp =  abs(varOriCell(1,k,m)-varOriCellNearCell(1,ikesi,ieta,m))
                Vtemp_y(m) = max(Vtemp_y(m),temp)
        
                do l = 2,nsp
                    temp =  abs(varOriCell(l,k,m)-varOriCell(l-1,k,m))
                    Vtemp_y(m) = max(Vtemp_y(m),temp)
                end do 
                !y���򣬿������3�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
                if(sideth(3)==1)then    
                    ikesi = 1;          ieta = k;
                elseif(sideth(3)==2)then
                    ikesi = k;    ieta = nsp
                elseif(sideth(3)==3)then
                    ikesi = nsp;   ieta = nsp+1-k;
                elseif(sideth(3)==4)then
                    ikesi = nsp+1-k;     ieta = 1
                end if  
                temp =  abs(varOriCell(nsp,k,m)-varOriCellNearCell(3,ikesi,ieta,m))
                Vtemp_y(m) = max(Vtemp_y(m),temp)
                
                !������ⵥԪ
                if(Vtemp_x(m)>theta*Max_TV(m))then
                    cellset(i).Beta_line(k,1) = 1
                    cellset(i).Beta = 1
                    if(detect_type == ByCell )then
                        cellset(i).Beta_line = 1
                        exit   
                    end if 
                end if
                if(Vtemp_y(m)>theta*Max_TV(m))then
                    cellset(i).Beta_line(k,2) = 1
                    cellset(i).Beta = 1
                    if(detect_type == ByCell )then
                        cellset(i).Beta_line = 1
                        exit   
                    end if
                end if
            end do     
        end do        
    end do
    

end subroutine indicator_MV