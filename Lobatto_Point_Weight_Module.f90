module Lobatto_P_W
    
    !-----------------------------------------------------------------------------
    !
    !   ������Gauss-Lobatto���ָ�˹�㼰Ȩ�ص����ģ��
    !   ������������������Lobatto���Ȩֵ����ʵ������߽�CPR �����̵档��Ҫ������������Զ����ɺ�������δ��ӡ� 
    !   �ο���https://keisan.casio.com/exec/system/1280801905
    !   ���ߣ�gqShi 
    !   ��ʷ������������δʵ���������� //2021.12.14
    !
    !-----------------------------------------------------------------------------    
    
    use real_precision
    use parameter_setting,only:nsp
     
    implicit none 
    
    integer :: n
    real(prec_dd),allocatable :: Lobatto_Point(:),Lobatto_Weight(:)
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
        
        n = nsp+1                                             !��������˹��ĸ���
        allocate(Lobatto_Point(n),Lobatto_Weight(n))
    end subroutine get_pre
    
    function f_Legendre_n(x,p)                                   !���� p-1 �׾��ȵĺ���f_Legendre_n(x)
        implicit none
        integer::i,p
        real(prec_dd) :: a(p),x,f_Legendre_n                 !a(p)����p�����õ¶���ʽ
        a(1) = x                                    !1�����õ¶���ʽ
        a(2) = 1.5_prec*(x**2)-0.5_prec                 !2�����õ¶���ʽ
        if(p>=3)then
            do i = 3,p                                  !���õ��ƹ�ϵ����p�����õ¶���ʽ
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
            end do            
        end if
    
        f_Legendre_n = a(p)                                    !���ɵ� p �����õ¶���ʽ 
    end function

    function f_Legendre_derivative(x,p)!����n�����õ¶���ʽ�ĵ������ʽ
        implicit none
        integer :: i,p
        real(prec_dd) :: a(p), x, f_Legendre_derivative
        a(1) = x
        a(2) = 1.5_prec*x**2-0.5_prec

        if(p>=3)then
            do i = 3,p                                  !���õ��ƹ�ϵ����p�����õ¶���ʽ
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
            end do            
        end if
        f_Legendre_derivative = p*a(p-1)/(1-x**2)-p*x*a(p)/(1-x**2)         
    end function

    function f_dichotomy(a,b,p)                   !���ַ���⺯���Ľ�
        implicit none
        real(prec_dd)::a,b,c,f_dichotomy                          
        integer :: p
        
        do while(.true.)
            c = (a+b)/2.0_prec
            if(f_Legendre_derivative(c,p-1)*f_Legendre_derivative(a,p-1) < 0)then
                b = c
            elseif(f_Legendre_derivative(c,p-1)*f_Legendre_derivative(a,p-1) > 0)then
                a = c
            else
                exit
            endif
            if((b-a) < eps)exit
            
        enddo
        f_dichotomy = c         
    end function
    
    subroutine over_post
        implicit none
        deallocate(Lobatto_Point,Lobatto_Weight)
    end subroutine over_post
        
    subroutine get_Lobatto_Point_Weight(Lobatto_Point,Lobatto_Weight)
        implicit none
        real(prec_dd) :: m,delta,Lobatto_Point(:),Lobatto_Weight(:)    !�Ѿ�����            
        integer :: i,j

        call get_pre

        j = 1                                                                                 
        Lobatto_Point(j) = -1.0_prec
        Lobatto_Weight(j) =  2.0_prec/(n*(n-1)*(f_Legendre_n(Lobatto_Point(j),n-1))**2)  
        m = -1.0001_prec                                                                                        
        do i = 1,20000                                                                                          
            if(f_Legendre_derivative(m,n-1)*f_Legendre_derivative(m+0.0001_prec ,n-1)<0)then               
                j = j+1                                                                                 
                Lobatto_Point(j) = f_dichotomy(m,m+0.0001_prec ,n)                                         
                                        
                Lobatto_Weight(j) =  2.0_prec/(n*(n-1)*(f_Legendre_n(Lobatto_Point(j),n-1))**2)            
            endif
             m = m+0.0001_prec                                 
        enddo
        j = j+1 
        Lobatto_Point(j) = 1.0_prec
        Lobatto_Weight(j) =  2.0_prec/(n*(n-1)*(f_Legendre_n(Lobatto_Point(j),n-1))**2)
        !write(*,*) ' '
        !
        !do j = 1,n         
        !    write(*,*)  'Lobatto��',Lobatto_Point(j)
        !end do
        !write(*,*) ' '
        !do j = 1,n         
        !    write(*,*)  'LobattoȨ',Lobatto_Weight(j)
        !end do
        !

        call over_post
    end subroutine
end module Lobatto_P_W 