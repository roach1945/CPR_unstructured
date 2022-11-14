program main

    !-----------------------------------------------------------------------------
    !
    !   方法：非结构四边形网格CPR方法的子单元限制技术
    !   目的：测试精度，捕捉间断
    !   描述：网格生成，网格读取，网格变换，CPR计算，精度测试，子单元限制
    !   作者：gqShi
    !   版本：v1.0 
    !   历史：重新梳理注释 //2021.12.08
    !
    !-----------------------------------------------------------------------------
    
    
    use preParameter
    use parameter_setting
    use global_var
    use Time_module
    use type_module
    implicit none
    
    ! 读取预设参数
    
    call read_preParameter      
    
    !-----------------------------------------------------------------------------
    !
    !   描述：
    !         计算过程：预处理-时间推进-求解准确解-后处理  //计算时间统计 
    !         功能实现：流场计算、收敛精度计算（自生成网格，网格尺寸成倍数减小）
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
    
    select case(compuPurpose)
    case(error_solve)
        
        ! 流场计算
        
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
        
        ! 误差收敛阶计算 10 20 30 40 ...
        
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
        
        ! 误差收敛阶计算 10 20 40 80 ...
        
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
