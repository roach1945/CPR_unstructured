program main

    !-----------------------------------------------------------------------------
    !
    !   �������ǽṹ�ı�������CPR�������ӵ�Ԫ���Ƽ���
    !   Ŀ�ģ����Ծ��ȣ���׽���
    !   �������������ɣ������ȡ������任��CPR���㣬���Ȳ��ԣ��ӵ�Ԫ����
    !   ���ߣ�gqShi
    !   �汾��v1.0 
    !   ��ʷ����������ע�� //2021.12.08
    !
    !-----------------------------------------------------------------------------
    
    
    use preParameter
    use parameter_setting
    use global_var
    use Time_module
    use type_module
    implicit none
    
    ! ��ȡԤ�����
    
    call read_preParameter      
    
    !-----------------------------------------------------------------------------
    !
    !   ������
    !         ������̣�Ԥ����-ʱ���ƽ�-���׼ȷ��-����  //����ʱ��ͳ�� 
    !         ����ʵ�֣��������㡢�������ȼ��㣨��������������ߴ�ɱ�����С��
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
    
    select case(compuPurpose)
    case(error_solve)
        
        ! ��������
        
        call cpu_time(TimeStart)
        call pre_processing         
        call cpu_time(TimeFrame1)
        call time_advancing     
        call cpu_time(TimeFrame2)       
        call exact_solution(T)     
        call cpu_time(TimeFrame3)
        call post_processing
        call cpu_time(TimeEnd)
        call Time_com
        
    case(order_solve)
        
        ! ��������׼��� 10 20 30 40 ...
        
        nsdx = 10
        nsdy = nsdx
        do while(nsdx<=50)
            call pre_processing     
            call time_advancing    
            call exact_solution(T)  
            call post_processing
            call solve_error_order                     
            
            nsdx = nsdx+10
            nsdy = nsdx            
        end do
        
    case(order_solve2)
        
        ! ��������׼��� 10 20 40 80 ...
        
        nsdx = 10
        nsdy = nsdx
        do while(nsdx<=80)
            call pre_processing     
            call time_advancing     
            call exact_solution(T)  
            call post_processing    
            call solve_error_order2
            
            nsdx = nsdx*2
            nsdy = nsdx           
        end do
        
    case default
    
        write(*,*) 'Non computed!'
        
    end select
    
    
    close(1000)
end program main
