!
! Fortran version of Little Boxes Simple 
! Mostly for making the plots much quicker for my dissertation 
!
!

program main 

    implicit none

    ! Declare general variables and parameters 

    real (kind=8), parameter :: start_time = 0.0d0
    complex (kind=8), dimension(2) :: coeffs 
    integer, parameter :: end_time = 20d0 
    integer (kind=8), parameter :: time_steps = end_time * 5000d0 
    real (kind=8), parameter :: dt = real(end_time) / real(time_steps)
    real (kind=8), dimension(time_steps) :: time_list 
    integer (kind=8), parameter :: num_of_simulations = 1000d0 
    complex (kind=8), parameter :: Omega = 2.3d0 
    real (kind=8), parameter :: pi = 3.14159265358979323846d0 
    real (kind=8) :: total, rand_num 
    integer :: beginning, end, rate, t, sim, index 
    complex (kind=8) :: g_0, e_0, g_1, e_1, g_2, e_2
    complex (kind=8) :: g_0_new, e_0_new, g_1_new, e_1_new, g_2_new, e_2_new
    real (kind=8), dimension(time_steps) :: rand_list 

    ! sigma z, lowering, raising 
    real (kind=8), dimension(time_steps) :: avg_sigma_z_list, avg_sigma_L_list, avg_sigma_R_list, g2

    ! Start running code 

    print *, dt 
    
    ! Initialise arrays 
    avg_sigma_L_list = 0.0d0; avg_sigma_R_list = 0.0d0; avg_sigma_z_list = 0.0d0; g2 = 0.0d0 

    ! Initialise time list array 
    call linspace(start=start_time, end=end_time, time_list=time_list)

    ! Program execution time tracking 
    call system_clock(beginning, rate)

    do sim = 1, num_of_simulations

        ! INITIALISE ARRAYS HERE AND CLEAR ARRAYS HERE 
        coeffs = 0.0d0; coeffs(1) = 1.0d0 
        g_0_new = 0.0d0; e_0_new = 0.0d0; g_1_new = 0.0d0; e_1_new = 0.0d0

        call random_number(rand_list)

        do t = 1, size(time_list)

            g_0_new = coeffs(1) - ((Omega/2)*coeffs(2)*dt)
            e_0_new = (coeffs(2)*(1-dt/2)) + (Omega/2)*coeffs(1)*dt
            g_1_new = cmplx(0,-1) * sqrt(dt/2) * coeffs(2)
            e_1_new = cmplx(0,-1) * sqrt(dt/2) * coeffs(2)

            ! if (t == 1) then 
            !     print *, dt 
            !     print *, g_0_new 
            !     print *, e_0_new 
            !     print *, coeffs(1)
            !     print *, coeffs(2)
            ! end if 

            total = sqrt(modulo_func(e_0_new)**2 + &
                    modulo_func(g_0_new)**2 + modulo_func(e_1_new)**2 &
                    + modulo_func(g_1_new)**2)

            ! I THINK IT SHOULD BE NOT SQUARED...

            g_0_new = g_0_new / total 
            e_0_new = e_0_new / total 
            g_1_new = g_1_new / total 
            e_1_new = e_1_new / total 

            rand_num = rand_list(t)

            if (rand_num <= 2 * modulo_func(g_1_new)**2) then 

                ! Update the sigma z list 
                avg_sigma_z_list(t) = avg_sigma_z_list(t) - 1.0d0 

                g2(t) = g2(t) + (g_0_new * e_0_new)
                
                ! update coeffs list 
                coeffs = 0.0d0
                coeffs(1) = 1.0d0 

            else

                total = sqrt(modulo_func(e_0_new)**2 + &
                    modulo_func(g_0_new)**2)

                g_0_new = g_0_new / total 
                e_0_new = e_0_new / total 

                ! Update the sigma z, L, R lists 
                avg_sigma_z_list(t) = avg_sigma_z_list(t) + modulo_func(e_0_new)**2 - modulo_func(g_0_new)**2 
                avg_sigma_L_list(t) = avg_sigma_L_list(t) + complex_multiply(conjg(g_0_new),e_0_new)
                avg_sigma_R_list(t) = avg_sigma_R_list(t) + complex_multiply(conjg(e_0_new),g_0_new)

                g2(t) = g2(t) + (g_0_new * e_0_new)

                ! Update coeffs list 
                coeffs = 0.0d0 
                coeffs(1) = g_0_new 
                coeffs(2) = e_0_new  

            end if 

        end do 

    end do 


    ! Normalise arrays 
    avg_sigma_z_list = avg_sigma_z_list / num_of_simulations
    avg_sigma_L_list = avg_sigma_L_list / num_of_simulations
    avg_sigma_R_list = avg_sigma_R_list / num_of_simulations
    g2 = g2 / num_of_simulations

    !!! Write out final results to a txt file 
    open(1, file="sigma_z.txt", status="replace")
    open(2, file="sigma_L.txt", status="replace")
    open(3, file="sigma_R.txt", status="replace")   
    open(4, file="g2.txt", status="replace")

    do index = 1,size(time_list)
        write(1,*) time_list(index), avg_sigma_z_list(index)
        write(2,*) time_list(index), avg_sigma_L_list(index)
        write(3,*) time_list(index), avg_sigma_R_list(index)
        write(4,*) time_list(index), g2(index)
    end do 

    close(1); close(2); close(3)

    call system_clock(end)

    print *, "All simulations completed. Execution time: ", real(end - beginning) / real(rate), " seconds."


    ! -------------------------------------------------------------------------------------------
    ! 
    !   Functions and Subroutines 
    !
    !-------------------------------------------------------------------------------------------

    contains 

    subroutine linspace(start, end, time_list) 

        real (kind=8), intent(in) :: start
        integer, intent(in) :: end 
        real (kind=8), intent(out) :: time_list(:)
        real (kind=8) :: range
        integer :: n, i

        n = size(time_list)
        range = end - start
        
        do i = 1,n
            time_list(i) = start + (range * (i - 1) / (n - 1))
        end do

    end subroutine

    function modulo_func(z) result(c)

        ! Takes in a complex number z and returns its modulus

        implicit none

        ! Declare var types
        real (kind=8) :: a, b, c
        complex (kind=8), intent(in) :: z

        a = real(z)
        b = aimag(z)

        c = sqrt(a**2 + b**2)

    end function 

    function complex_multiply(a,b) result(c)

        ! a is the complex conjugated variable 

        implicit none 

        real (kind=8) :: c
        complex (kind=8), intent(in) :: a,b 

        c = real(a)*real(b) + aimag(a)*aimag(b)

    end function

    
end program main 