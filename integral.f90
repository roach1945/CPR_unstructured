
function Gauss_integral(xl,xr,func,order)
    use real_precision
    use global_var,only:GPs,Gcoe
    real(prec),external :: func
    real(prec) :: xl,xr,Gauss_integral
    integer :: order,i,j 

    lamta1 = (b - a)*0.5_prec
    lamta2 = (b + a)*0.5_prec
    
    Gauss_integral = 0.0_prec

    do i=1,order
        temp = lamta2 + lamta1 * GPs(i)
        u_trans = lamta1 * func(temp)
        
        Gauss_integral = Gauss_integral + GCoe(i) * u_trans!积分公式
    end do
    return
end function Gauss_integral
function Gauss_integral_SPs(varFunc)
    !计算标准单元积分，利用已有的Gauss点，Gauss权
    !若解点不为Gauss点，需要先插值到Gauss点上。暂记
    use parameter_setting
    use global_var
    implicit none
    real(prec),dimension(nsp) :: varFunc
    real(prec) :: temp,temp2,Gauss_integral_SPs
    integer i,j
    temp = 0.0_prec
    do j = 1,nsp
        temp = temp + Gcoe(j)*varFunc(j)
    end do    
    Gauss_integral_SPs = temp
    return

end function Gauss_integral_SPs
function Gauss_double_integral(varFunc)
    !计算标准单元积分，利用已有的Gauss点，Gauss权
    use parameter_setting
    use global_var
    implicit none
    real(prec),dimension(nsp,nsp) :: varFunc
    real(prec) :: temp,temp2,Gauss_double_integral
    integer i,j
    temp2 = 0.0_prec
    do i = 1,nsp
        temp = 0.0_prec
        do j = 1,nsp
            temp = temp + Gcoe(j)*varFunc(i,j)
        end do
        temp2 = temp2 + Gcoe(i)*temp
    end do
    
    Gauss_double_integral = temp2
    return
end function Gauss_double_integral

subroutine Gauss_DInt(varFunc,Gauss_DInt_result)
    !计算标准单元积分，利用已有的Gauss点，Gauss权
    use parameter_setting
    use global_var
    implicit none
    real(prec),dimension(nsp,nsp) :: varFunc
    real(prec) :: temp,temp2,Gauss_DInt_result
    integer i,j
    Gauss_DInt_result = 0.0_prec
    temp2 = 0.0_prec
    do i = 1,nsp
        temp = 0.0_prec
        do j = 1,nsp
            temp = temp + Gcoe(j)*varFunc(i,j)
        end do
        temp2 = temp2 + Gcoe(i)*temp
    end do
    
    Gauss_DInt_result = temp2
end subroutine Gauss_DInt