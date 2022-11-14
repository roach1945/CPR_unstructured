! ���޸�
module Correction_Func  
     
    !-----------------------------------------------------------------------------
    !
    !   �������������ϵ��ģ��
    !   ���������ݲ�ͬ������������������ǰ�����Ӧ������ϵ������gDG��g2��
    !   ���ߣ�gqShi 
    !   ��ʷ������ //2021.12.15
    !
    !-----------------------------------------------------------------------------
    
    use real_precision
    use parameter_setting,only:nsp

    implicit none 
    integer :: n                        
    real(prec_dd),allocatable :: gl_coefficient(:),gr_coefficient(:),SPs_dd(:)  
    real(ddouble_prec) :: eps
    
    contains
    
    subroutine get_pre
        implicit none
        if(prec_dd == ddouble_prec)then
            eps = 1E-34
        elseif(prec_dd == double_prec )then
            eps =  1E-15
        elseif(prec_dd == single_prec )then
            eps = 1E-8   
        end if
        
        n = nsp                                             
        allocate(gl_coefficient(n),gr_coefficient(n),SPs_dd(n))
    end subroutine get_pre
    
    subroutine over_post
        implicit none
        deallocate(gl_coefficient,gr_coefficient,SPs_dd)
    end subroutine over_post
        
    function f_Legendre_derivative(x,p)!����n�����õ¶���ʽ�ĵ������ʽ
        implicit none
        integer :: i,p
        real(prec_dd) :: a(p), x, f_Legendre_derivative
        a(1) = x
        a(2) = 1.5_prec_dd*x**2-0.5_prec_dd

        if(p>=3)then
            do i = 3,p
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
                !write(*,*)a(i)
            enddo            
        end if

        f_Legendre_derivative = p*a(p-1)/(1.0_prec_dd-x**2)-p*x*a(p)/(1.0_prec_dd-x**2)     
        
    end function f_Legendre_derivative
     
    subroutine get_gDG_collocation(gl_coefficient,gr_coefficient)
        use global_var,only: SPs
        implicit none
        real(prec_dd) :: gl_coefficient(:),gr_coefficient(:)               
        integer :: i
        
        call get_pre
        
        do i = 1,n         
            SPs_dd(i) =  SPs(i)
        end do
                                                                                      
        do i = 1,n   
            gl_coefficient(i) = (-1.0_prec_dd)**n*0.5_prec_dd*(f_Legendre_derivative(SPs_dd(i),n)-f_Legendre_derivative(SPs_dd(i),n-1))
            gr_coefficient(i) = 0.5_prec_dd*(f_Legendre_derivative(SPs_dd(i),n)+f_Legendre_derivative(SPs_dd(i),n-1))        
        enddo
        
        call over_post
     
    end subroutine get_gDG_collocation
    
    
    
    !subroutine get_sp_collocation(gl_coefficient,gr_coefficient)
    !    ! ����������������ķ���.Ŀǰδ�ӣ������ӿ�
    !
    !    use parameter_setting,only:corFunv_type,corFunv_gDG !δ���壬ʵ���뷨��
    !    implicit none
    !    
    !    real(prec_dd) :: gl_coefficient(:),gr_coefficient(:)  
    !    
    !    select case(corFunv_type)
    !    case(corFunv_gDG)
    !        call get_gDG_collocation(gl_coefficient,gr_coefficient)
    !    case(corFunv_g2)
    !        
    !    case default
    !        
    !    end seclect
    !               
    !end subroutine get_sp_collocation
    
end module Correction_Func

    
    
    
    
! �ɰ汾�������ϵ���ķ���������ʵ�����������ϵ�����ģ�飬�ʿ����á�    
function gl_sub_kesi(kesi)
    !!gl = (-1)**K/2*(P_k-P_k-1),Pk is Legendre polynomial,��˽�������
    !!K=1,gl = -1/2*(x-1),gl_sub_kesi = -1/2
    !!K=2,gl = 1/4*(3*x**2-2*x-1),gl_sub_kesi = 1/4*(6*x-2)
    !!K=3,gl = -1/4*(5*x**3-3*x**2-3*x+1),gl_sub_kesi = -1/4*(15*x**2-6*x-3)
    !!K=4,gl = 1/16*(35*x**4-20*x**3-30*x**2+12*x+3),gl_sub_kesi = 1/16*(140*x**3-60*x**2-60*x+12)
    !!K=5,gl = -1/16*(63*x**5-35*x**4-70*x**3+30*x**2+15*x-3),gl_sub_kesi = -1/16*(315*x**4-140*x**3-210*x**2+60*x+15)
    use real_precision
    use parameter_setting ,only:nsp
    implicit none
    real(prec) :: gl_sub_kesi,kesi
    if (nsp .EQ. 1)then
        gl_sub_kesi = -0.5_prec
    elseif (nsp .EQ. 2)then
        gl_sub_kesi = 0.25*(6.0_prec*kesi - 2.0_prec)
    elseif (nsp .EQ. 3)then
        gl_sub_kesi = -0.25_prec*(15.0_prec*kesi**2.0_prec-6.0_prec*kesi-3.0_prec)
    elseif(nsp .EQ. 4)then
        gl_sub_kesi = (140.0_prec*kesi**3-60.0_prec*kesi**2-60.0_prec*kesi+12.0_prec)/16.0_prec
    elseif(nsp .EQ. 5)then
        gl_sub_kesi = (315.0_prec*kesi**4-140.0_prec*kesi**3-210.0_prec*kesi**2+60.0_prec*kesi+15.0_prec)*(-1.0_prec)/16.0_prec
    endif
    return
end function gl_sub_kesi

function gr_sub_kesi(kesi)
    !!gr = 1/2*(P_k+P_k-1),�Ҷ˽�������
    !!K=1,gr = 1/2*(x+1),gr_sub_kesi = 1/2
    !!K=2,gr = 1/4*(3*x**2+2*x-1),gl_sub_kesi = 1/4*(6*x+2)
    !!K=3,gr = 1/4*(5*x**3+3*x**2-3*x-1),gr_sub_kesi = 1/4*(15*x**2+6*x-3)
    !!K=4,gr = 1/16*(35*x**4+20*x**3-30*x**2-12*x+3),gr_sub_kesi = 1/16*(140*x**3+60*x**2-60*x-12)
    !!K=5,gr = 1/16*(63*x**5+35*x**4-70*x**3-30*x**2+15*x+3),gr_sub_kesi = 1/16*(315*x**4+140*x**3-210*x**2-60*x+15)
    use real_precision
    use parameter_setting ,only:nsp
    implicit none
    real(prec) :: gr_sub_kesi,kesi

    if (nsp .EQ. 1)then
        gr_sub_kesi = 0.5_prec
    elseif (nsp .EQ. 2)then
        gr_sub_kesi = 0.25*(6.0_prec*kesi + 2.0_prec)    
    elseif (nsp .EQ. 3)then
        gr_sub_kesi = 0.25_prec*(15.0_prec*kesi**2.0_prec+6.0_prec*kesi-3.0_prec)
    elseif(nsp .EQ. 4)then
        gr_sub_kesi = (140.0_prec*kesi**3+60.0_prec*kesi**2-60.0_prec*kesi-12.0_prec)/16.0_prec
    elseif(nsp .EQ. 5)then
        gr_sub_kesi = (315.0_prec*kesi**4+140.0_prec*kesi**3-210.0_prec*kesi**2-60.0_prec*kesi+15.0_prec)/16.0_prec
    endif
    return
end function gr_sub_kesi