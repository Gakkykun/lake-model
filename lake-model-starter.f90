program lake
    implicit none

    real, allocatable :: cChla(:, :), cTOC(:, :), cTP(:, :), cTN(:, :), cDO(:, :), cSS(:, :), t(:, :)
    integer :: q, i, j, N, day, month, count_day
    real :: dt, tt, L, dz
    real :: dAdy, B, alpha, diff, rho, C_heat
    real :: Gy_max, Ry_max, kc_20, kb_20, theta, kp, Is, w0, &
            gamma_cy, gamma_py, gamma_ny, beta_y, gamma_oy, gamma_oc, yeta
    real :: elevation, Q_in, Q_out, V, phi_n, T_in, Chla_in, TOC_in, TP_in, TN_in, &
            DO_in, SS_in, u_in, u_out, Vw, month_c
    real :: A, phi_y, fi, fp, Gy, Ry, kc, kb
    real :: s1, s2, v1Chla, v2Chla, v1TOC, v2TOC, v1TP, v2TP, v1TN, v2TN, &
            v1DO, v2DO, v1SS, v2SS, v1T, v2T
    real, dimension(61, 5) :: MonitoringDataDaily
    real, dimension(12, 7) :: MonitoringDataMonth

    dt = 0.2                ! Time step [d]
    tt = 61.2               ! Total time [d]
    q = int(tt/dt)          ! Total time steps
    L = 80.0                ! Total depth [m]
    N = 40                  ! Total space steps
    dz = L/real(N)          ! Space step: 0.2 [m]


    ! Parameters
    ! advection and dispersion
    B = 800.0                   ! Width of each layer [m]
    dAdy = 62.5 * B             ! variation of surface area [m^2]   
    alpha = 9.0*10.0**(-4)      ! (Thermal )Eddy diffusion coefficient [m^2/d]
    diff = 9.0*10.0**(-5)       ! Molecular (heat) diffusion coefficient [m^2/d]
    rho = 1.0                   ! Water density [t/m^3]
    C_heat = 1.0                ! Heat capacity [kcal/kg/C]  

    s1 = dt / dz
    s2 = (alpha + diff) * dt / dz**2

    ! process rates in mass balance equation
    Gy_max = 2.0                ! Maximum specific growth rate of algae at 20 degree of C[1/d]
    Ry_max = 0.09               ! Endogenous respiration rate of algae at 20 degree of C [1/d]
    kc_20 = 0.01                ! Decomposition rate constant of TOC at 20 degree of C [1/d]
    kb_20 = 2.5                 ! Oxygen utilization rate of sediment at 20 degree of C [1/d] 
    theta = 1.07                ! Temperature coefficient [-]
    kp = 0.026                  ! Saturation coefficient [mg/L]
    Is = 300                    ! Optimum light intensity [kcal/m^2/d]
    w0 = 0.4                    ! Settling velocity [m/d]
    gamma_cy = 4.8*10.0**(-2)   ! Rate of biological release of TOC from algae[gTOC/mgChla]
    gamma_py = 1.25             ! Coefficient of biological utilization of TP for algae [mgP/mgChla]
    gamma_ny = 17.0             ! Coefficient of biological utilization of TN for algae [gTN/mgChla]
    beta_y = 0.6                ! Yield coefficient of TN and TP[-]
    gamma_oy = 0.13             ! Coefficient of biological production of O2 for algae [gO2/mgChla]
    gamma_oc = 2.7              ! Oxygen demand exerted by TOC decomposition [gO2/gTOC]
    yeta = 0.4                  ! Lambert's law


    ! 1. Allocation of array
    
    
    ! 2. Initial and boundary conditions
    

    ! 3. Import monitoring data

    ! Loop calculation
    do j = 1, q-1
        ! 4. Monitoring data - daily average
                        
            ! Horizontal and vertical velocities
            u_in = Q_in/(B*L)        ! Inlet horizontal velocity [m/d]
            u_out = Q_out/(B*L)      ! Outlet horizontal velocity [m/d] 
            V = Q_out/A              ! Vertical velocity [m/d]
            
            
            ! 4.1 Monitoring data - montyly average
          
        

    do i = 1, N

        ! 5. Calculations of parameters included in governing equations
        

        ! 6. Chl.a
         
        
        if (cChla(i, j+1) < 0 ) then
            cChla(i, j+1) = 0
        endif
        
        ! 7. TOC

        if (cTOC(i, j+1) < 0 ) then
            cTOC(i, j+1) = 0
        endif
        
        ! 8. TP
        
        if (cTP(i, j+1) < 0 ) then
            cTP(i, j+1) = 0
        endif
        
        ! 9. TN
        
        if (cTN(i, j+1) < 0 ) then
            cTN(i, j+1) = 0
        endif
        
        ! 10. DO
        
        if (cDO(i, j+1) < 0 ) then
            cDO(i, j+1) = 0
        endif

        ! 11. SS

        if (cSS(i, j+1) < 0 ) then
            cSS(i, j+1) = 0
        endif
        
        ! 12. Temperature
        
        if (t(i, j+1) < 0 ) then
            t(i, j+1) = 0
        endif

    end do ! for i (position)

    open(10,file='output_chla.csv')
    open(11,file='output_toc.csv')
    open(12,file='output_tp.csv')
    open(13,file='output_tn.csv')
    open(14,file='output_do.csv')
    open(15,file='output_ss.csv')
    open(16,file='output_t.csv')

    ! 13. Writing results into csv files

    end do ! for j (time)

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    contains    ! 14. Two functions: light_intensity and aeration
    
  


end program lake
