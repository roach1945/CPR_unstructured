module Gauss_P_W

    !-----------------------------------------------------------------------------
    !
    !   ��������˹�����õ»��ָ�˹�㼰Ȩ�ص����ģ��
    !   ������������������Gauss���Ȩֵ����ʵ������߽�CPR �����̵档��Ҫ������������Զ����ɺ�������δ��ӡ� 
    !   �ο���http://bbs.fcode.cn/thread-219-1-1.html 
    !   ���ߣ�gqShi 
    !   ��ʷ������������δʵ���������� //2021.12.14
    !
    !-----------------------------------------------------------------------------    
    
    use real_precision
    use parameter_setting,only:nsp

    
    implicit none 
    integer :: n                        
    real(prec_dd),allocatable :: Gauss_Point(:),Gauss_Weight(:)
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
        
        n = nsp                                             !��������˹��ĸ���
        allocate(Gauss_Point(n),Gauss_Weight(n))
    end subroutine get_pre
    
    subroutine over_post
        implicit none
        deallocate(Gauss_Point,Gauss_Weight)
    end subroutine over_post
    
    function f_Legendre_n(x)                                !���� n-1 �׾��ȵĺ���f_Legendre_n(x)
        implicit none
        integer::i
        real(prec_dd) :: a(n),x,f_Legendre_n                !a(n)����n�����õ¶���ʽ
        a(1) = x                                            !1�����õ¶���ʽ
        a(2) = 1.5_prec*(x**2)-0.5_prec                     !2�����õ¶���ʽ

        if(n>=3)then
            do i = 3,n
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
            enddo            
        end if
        f_Legendre_n = a(n)                                 !���ɵ� n �����õ¶���ʽ 
    end function

    function f_Legendre_nMinusOne(x) 
        implicit none
        integer::i
        real(prec_dd)::a(n), x, f_Legendre_nMinusOne
        a(1) = x
        a(2) = 1.5_prec*x**2-0.5_prec
        do i = 3,n-1
            a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
        enddo
        f_Legendre_nMinusOne = a(n-1)                                   !���ɵģ�n-1)�����õ¶���ʽ 
    end function

    function f_Legendre_derivative(x)
        implicit none
        integer :: i
        real(prec_dd) :: a(n), x, f_Legendre_derivative
        a(1) = x
        a(2) = 1.5_prec*x**2-0.5_prec
        if(n>=3)then
            do i = 3,n
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
                !write(*,*)a(i)
            enddo            
        end if
        f_Legendre_derivative = n*a(n-1)/(1-x**2)-n*x*a(n)/(1-x**2)     !����n�����õ¶���ʽ�ĵ������ʽ.������ƹ�ʽ
    end function

    function f_dichotomy(a,b)                                           !���ַ���⺯���Ľ�
        implicit none
        real(prec_dd)::a,b,c,f_dichotomy                                !a,b�Ǵ��ݽ����Ļ��ֺõ���һ������ڵ�����
        do while(.true.)
            c = (a+b)/2.0_prec
            if(f_Legendre_n(c)*f_Legendre_n(a) < 0)then
                b = c
            else
                a = c
            endif
            if((b-a) < eps)exit
            
        enddo
        f_dichotomy = c                                                     !f_dichotomy�������ö��ַ���õĽ�
    end function

    subroutine get_Gauss_Point_Weight(Gauss_Point,Gauss_Weight)
        implicit none
        real(prec_dd) :: m,Gauss_Point(n),Gauss_Weight(n)                   !��������,��Сn��module��ʼ������
        integer :: i,j
        call get_pre

        j = 0                                                              !��ֵ����ѭ�������ĳ�ֵ 
        m = -1.00001                                                       !���ü�����[-1��1] �����ޣ�������-1 
        do i = 1,20000                                                     !���ѭ������Ӧ�����ɲ���0.00001�� ��,���㷽����200000=2/0.000001 
            if(f_Legendre_n(m)*f_Legendre_n(m+0.0001)<0)then               !�����޴���ʼ�������ۼӣ�!�ɲ���0.00001˵��������10^5����
                j = j+1                                                    !��¼���ǵڼ�����
                Gauss_Point(j) = f_dichotomy(m,m+0.0001)                   !���ö��ַ��������ڷֺõ�һС������⣬����洢��Gauss_Point��j)                                        
                Gauss_Weight(j) = 2.0_prec/(n*f_Legendre_nMinusOne(Gauss_Point(j))*f_Legendre_derivative(Gauss_Point(j)))     !��˹���Ȩ��
            endif
             m = m+0.0001                                
        enddo
        !do j = 1,n         
        !    write(*,*)  '��˹��',Gauss_Point(j)
        !end do
        !write(*,*) ' '
        !do j = 1,n         
        !    write(*,*)  '��˹Ȩ',Gauss_Weight(j)
        !end do
        call over_post

    end subroutine

end module Gauss_P_W