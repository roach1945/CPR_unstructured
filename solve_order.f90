subroutine Norm
    !计算范数误差
    use real_precision
    use type_module
    use global_var
    implicit none
    integer i,j,k,l
    real(prec) :: err,errM(nsp,nsp)
    real(prec),external :: Gauss_double_integral
    err = 0.0_prec
    Norm_L1 = 0.0_prec
    Norm_L2 = 0.0_prec
    Norm_Linf = 0.0_prec
    do i = 1,ncells
        do j = 1,nsp 
            do k =1,nsp          
                err = abs(cellset(i).spvalue_ori(j,k,1)-cellset(i).spvalue_ori_exa(j,k,1))
                Norm_L1 = Norm_L1 + err
                Norm_L2 = Norm_L2 + err**2
                Norm_Linf = max(Norm_Linf,err)
            end do
        end do   
    end do
    Norm_L1 = Norm_L1/(ncells*nsp*nsp) 
    Norm_L2 = sqrt(Norm_L2/(ncells*nsp*nsp))  
    
    !write(*,*)Norm_L1!/(ncells*nsp*nsp) 
    !Norm_L1 = 0.0_prec
    !do i = 1,ncells
    !    errM = abs(cellset(i).spvalue_ori(:,:,1)-cellset(i).spvalue_ori_exa(:,:,1))
    !    
    !    Norm_L1 = Norm_L1 + Gauss_double_integral(errM)*cellset(i).det_J(3,3)
    !    Norm_L2 = Norm_L2 + Gauss_double_integral(errM)**2
    !    Norm_Linf = max(Norm_Linf,MAXVAL(errM))
    !    write(*,*)cellset(i).det_J(3,3)
    !end do
    !Norm_L2 = sqrt(Norm_L2)  
    !write(*,*)Norm_L1
    !stop
 
    a_Norm_L1(counter) = Norm_L1
    a_Norm_L2(counter) = Norm_L2
    a_Norm_Linf(counter) = Norm_Linf
    counter = counter + 1
    
    
    write(1010,*) Norm_L1,Norm_L2,Norm_Linf
    !write(*,*) 'Norm_L1  ',Norm_L1
    !write(*,*) 'Norm_L2  ',Norm_L2
    !write(*,*) 'Norm_Linf',Norm_Linf
    write(*,*) Norm_L1
    write(*,*) Norm_L2
    write(*,*) Norm_Linf
end subroutine Norm


subroutine solve_error_order
    use global_var
    use parameter_setting
    implicit none
    integer i,j
    real(prec),dimension(10) :: order_L1,order_L2,order_Linf
    character(len=40) file_name
    character(len=16) char_T,char_nsd,char_scheme
    character(len=2) char_nsp
    character(len=2) char_ns

!!!!!!!!!求解范数误差 
    
    write(*,*) " "
    do i = 2,counter-1
        order_L1(i) = log(a_Norm_L1(i-1)/a_Norm_L1(i))/log(i/(i-1.0_prec))
        order_L2(i) = log(a_Norm_L2(i-1)/a_Norm_L2(i))/log(i/(i-1.0_prec))
        order_Linf(i) = log(a_Norm_Linf(i-1)/a_Norm_Linf(i))/log(i/(i-1.0_prec))
    end do
    write(*,*) "m       ","L1 error     ","L2 error     ","Linf error ","L1 order     ","L2 order   ","Linf order"
    do i = 1,counter-1
        write(*,"(1x,I3,3ES10.2,3F10.2 )") nsdx-10*(counter-i-1),a_Norm_L1(i),a_Norm_L2(i),a_Norm_Linf(i),order_L1(i),order_L2(i),order_Linf(i)
    end do
!!!!!!!!!Print error result to '*.dat' file    
    write(char_T,'(F16.4)') T !T写入为字符型char_T
    write(char_nsd,'(I4)') nsdx
    write(char_nsp,'(I2)') nsp
    
    select case(scheme_kind)
    case(scheme_cpr)
        char_scheme = 'cpr'
    case(scheme_two)
        char_scheme = 'cnnw2'
    case(scheme_hybrid)
        char_scheme = 'cpr-cnnw2'
    case default
        char_scheme = 'other'
    end select
    
    file_name = '.\error\'//trim(adjustl(char_scheme))//'_T='//trim(adjustl(char_T))//'__'//trim(adjustl(char_nsp))//'.txt'
    !inquire(file=file_name,exist=alive)
    open(unit=1001,file=file_name,position='append')
      
    i = counter-1
    if(i == 1)then
        write(1001,*) ' '
        write(1001,*) " m     ","L1 error       ","L2 error       ","Linf error       ","L1 order     ","L2 order   ","Linf order"
    end if
    
    write(1001,"(1x,I3,3ES10.2,3F10.2)") nsdx-10*(counter-i-1),a_Norm_L1(i),a_Norm_L2(i),a_Norm_Linf(i),order_L1(i),order_L2(i),order_Linf(i)

    close(1001)
end subroutine solve_error_order

subroutine solve_error_order2
    use global_var
    use parameter_setting
    implicit none
    integer i,j
    real(prec),dimension(10) :: order_L1,order_L2,order_Linf
    character(len=40) file_name
    character(len=16) char_T,char_nsd,char_scheme
    character(len=8) char_nsp
    character(len=8) char_ns

!!!!!!!!!求解范数误差 
    
    !write(*,*) " "
    !write(*,*) "   L1 error       ","   L2 error       ","Linf error"
    !do i = 1,counter-1
    !    write(*,*) a_Norm_L1(i),a_Norm_L2(i),a_Norm_Linf(i)
    !end do
    !write(*,*) "   L1 order       ","   L2 order       ","Linf order"
    do i = 2,counter-1
        order_L1(i) = log(a_Norm_L1(i-1)/a_Norm_L1(i))/log(2.0_prec)
        order_L2(i) = log(a_Norm_L2(i-1)/a_Norm_L2(i))/log(2.0_prec)
        order_Linf(i) = log(a_Norm_Linf(i-1)/a_Norm_Linf(i))/log(2.0_prec)
        !write(*,*) order_L1(i),order_L2(i),order_Linf(i)
    end do
    write(*,*) " m     ","L1 error       ","L2 error       ","Linf error       ","L1 order     ","L2 order   ","Linf order"
    do i = 1,counter-1
        write(*,"(1x,I3,3ES10.2,3F10.2)") 10*2**(i-1),a_Norm_L1(i),a_Norm_L2(i),a_Norm_Linf(i),order_L1(i),order_L2(i),order_Linf(i)
    end do
!!!!!!!!!Print error result to '*.dat' file    
    write(char_T,'(F6.4)') T !T写入为字符型char_T
    write(char_nsd,'(I4)') nsdx
    write(char_nsp,'(I2)') nsp
    
    select case(scheme_kind)
    case(scheme_cpr)
        char_scheme = 'cpr'
    case(scheme_two)
        char_scheme = 'cnnw2'
    case(scheme_hybrid)
        char_scheme = 'cpr-cnnw2'
    case default
        char_scheme = 'other'
    end select
    
    file_name = '.\error\'//trim(adjustl(char_scheme))//'_T='//trim(adjustl(char_T))//'__'//trim(adjustl(char_nsp))//'.txt'

    open(unit=1001,file=file_name,position='append')

    i = counter-1
    if(i == 1)then
        write(1001,*) ' '
        write(1001,*) " m     ","L1 error       ","L2 error       ","Linf error       ","L1 order     ","L2 order   ","Linf order"
    end if
    
    write(1001,"(1x,I3,3ES10.2,3F10.2)") 10*2**(i-1),a_Norm_L1(i),a_Norm_L2(i),a_Norm_Linf(i),order_L1(i),order_L2(i),order_Linf(i)
   
    close(1001)

end subroutine solve_error_order2