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


    ! Allocation of array
    allocate (cChla(0:N, 0:q))
    allocate (cTOC(0:N, 0:q))
    allocate (cTP(0:N, 0:q))
    allocate (cTN(0:N, 0:q))
    allocate (cDO(0:N, 0:q))
    allocate (cSS(0:N, 0:q))
    allocate (t(0:N, 0:q))
    
    ! Initial and boundary conditions
    cChla(:, 0) = 0.0
    cChla(0, :) = 0.0

    cTOC(:, 0) = 0.0
    cTOC(0, :) = 0.0

    cTP(:, 0) = 0.0
    cTP(0, :) = 0.0
    
    cTN(:, 0) = 0.0
    cTN(0, :) = 0.0
    
    cDO(:, 0) = 0.0
    cDO(0, :) = 0.0
    
    cSS(:, 0) = 0.0
    cSS(0, :) = 0.0
    
    t(:, 0) = 0.0
    t(0, :) = 0.0

    open(unit=99, file='MonitoringDataDaily.csv',status='old')
    read (99, '()') 
    open(unit=98, file='MonitoringDataMonth.csv',status='old')
    read (98, '()')

    day = 1
    month = 1
    count_day = int(1.0/dt) ! 5

    ! Loop calculation
    do j = 1, q-1
        ! Monitoring data - daily average
        if (mod(j, count_day) == 0) then
            day = day + 1
            
            read (99,*) MonitoringDataDaily(day, 1), MonitoringDataDaily(day, 2), &
                        MonitoringDataDaily(day, 3), MonitoringDataDaily(day, 4), &
                        MonitoringDataDaily(day, 5)
            
            Q_in = MonitoringDataDaily(day,1)              ! Inlet flow rate [m^3/d]
            Q_in = MonitoringDataDaily(day,1)              ! Inlet flow rate [m^3/d]
            Q_out = MonitoringDataDaily(day,2)             ! Outlet flow rate [m^3/d]         
            phi_n = MonitoringDataDaily(day,3)             ! Light intensity [kcal/m^3/d]
            Vw = MonitoringDataDaily(day,4)                ! Wind velocity [m/s]
            month_c = MonitoringDataDaily(day, 5)
                        
            ! Horizontal and vertical velocities
            u_in = Q_in/(B*L)        ! Inlet horizontal velocity [m/d]
            u_out = Q_out/(B*L)      ! Outlet horizontal velocity [m/d] 
            V = Q_out/A              ! Vertical velocity [m/d]
            

            ! Monitoring data - montyly average
            if (month == int(month_c)) then
                            
                read (98,*) MonitoringDataMonth(month, 1), MonitoringDataMonth(month, 2), &
                            MonitoringDataMonth(month, 3), MonitoringDataMonth(month, 4), &
                            MonitoringDataMonth(month, 5), MonitoringDataMonth(month, 6), &
                            MonitoringDataMonth(month, 7)
                            
                T_in = MonitoringDataMonth(month, 1)             ! Inlet water temperature 
                Chla_in = MonitoringDataMonth(month, 2)          ! Inlet chlorophyl concentration
                TOC_in = MonitoringDataMonth(month, 3)           ! Inlet TOC concentration
                TP_in = MonitoringDataMonth(month, 4)            ! Inlet TP concentration
                TN_in = MonitoringDataMonth(month, 5)            ! Inlet TN concentration
                DO_in = MonitoringDataMonth(month, 6)            ! Inlet DO concentration
                SS_in = MonitoringDataMonth(month, 7)            ! Inlet SS concentration
                
                month = month + 1
            end if
        end if

    do i = 1, N

        elevation = i*dz
        A = dAdy*elevation
        phi_y = light_intensity(elevation, phi_n)

        fi = phi_y/Is*exp(1-phi_y/Is)
        fp = cTP(i,j)/(kp+cTP(i,j))

        Gy = Gy_max*theta**(t(i,j)-20.0)*fi*fp
        Ry = Ry_max*theta**(t(i,j)-20.0)
        kc = kc_20*theta**(t(i,j)-20.0)
        kb = kb_20*theta**(t(i,j)-20.0)
        
        ! Chl.a
        v1Chla = s1/2*(cChla(i+1, j) - cChla(i-1, j))
        v2Chla = s2*(cChla(i+1, j) - 2*cChla(i,j) + cChla(i-1, j))

        cChla(i, j+1) = cChla(i, j) + dt*(Gy-Ry)*cChla(i,j) + W0*v1Chla + &
        dt*B/A*(u_in*Chla_in - u_out*cChla(i, j))- V*v1Chla + v2Chla

        if (i == N) then
            cChla(i, j+1) = cChla(i, j) + dt*(Gy-Ry)*cChla(i,j) + W0*v1Chla + &
            dt*B/A*(u_in*Chla_in - u_out*cChla(i, j)) + 2*s2*(cChla(i-1, j) - cChla(i, j))
        end if 
        
        if (cChla(i, j+1) < 0 ) then
            cChla(i, j+1) = 0
        endif
        
        ! TOC
        v1TOC = s1/2*(cTOC(i+1, j)- cTOC(i-1, j))
        v2TOC = s2*(cTOC(i+1, j) - 2*cTOC(i,j) + cTOC(i-1, j))

        cTOC(i, j+1) = cTOC(i, j) + dt*gamma_cy*Ry*cChla(i,j) - dt*kc*cTOC(i,j) + &
        dt*B/A*(u_in*TOC_in - u_out*cTOC(i, j))- V*v1toc + v2toc

        if (i == N) then
            cTOC(i, j+1) = cTOC(i, j) + dt*gamma_cy*Ry*cChla(i,j) - dt*kc*cTOC(i,j) + &
            dt*B/A*(u_in*TOC_in - u_out*cTOC(i, j)) + 2*s2*(cTOC(i-1, j) - cTOC(i, j))   
        end if

        if (cTOC(i, j+1) < 0 ) then
            cTOC(i, j+1) = 0
        endif
        
        ! TP
        v1TP = s1/2*(cTP(i+1, j)- cTP(i-1, j))
        v2Tp = s2*(cTP(i+1, j) - 2*cTP(i,j) + cTP(i-1, j))

        cTP(i, j+1) = cTP(i, j) + dt*gamma_py*(-Gy + beta_y*Ry)*cChla(i,j) + &
        dt*B/A*(u_in*TP_in - u_out*cTP(i, j))- V*v1Tp + v2Tp
        
        if (i == N) then
            cTP(N, j+1) = cTP(N-1, j+1) + dt*gamma_py*(-Gy + beta_y*Ry)*cChla(i,j) + &
            dt*B/A*(u_in*TP_in - u_out*cTP(i, j)) + 2*s2*(cTP(i-1,j) - cTP(i,j))
        end if
        
        if (cTP(i, j+1) < 0 ) then
            cTP(i, j+1) = 0
        endif
        
        ! TN
        v1TN = s1/2*(cTN(i+1, j)- cTN(i-1, j))
        v2TN = s2*(cTN(i+1, j) - 2*cTN(i,j) + cTN(i-1, j))

        cTN(i, j+1) = cTN(i, j) + dt*gamma_ny*(-Gy + beta_y*Ry)*cChla(i,j) + &
        dt*B/A*(u_in*TN_in - u_out*cTN(i, j))- V*v1Tn + v2Tn

        if (i == N) then
            cTN(i, j+1) = cTN(i, j) + + dt*gamma_ny*(-Gy + beta_y*Ry)*cChla(i,j) + &
            dt*B/A*(u_in*TN_in - u_out*cTN(i, j)) + 2*s2*(cTN(i-1, j) - cTN(i, j))   
        end if
        
        if (cTN(i, j+1) < 0 ) then
            cTN(i, j+1) = 0
        endif
        
        ! DO
        v1DO = s1/2*(cDO(i+1, j)- cDO(i-1, j))
        v2DO = s2*(cDO(i+1, j) - 2*cDO(i,j) + cDO(i-1, j))
        
        cDO(i, j+1) = cDO(i, j) + dt*gamma_oy*(Gy - Ry)*cChla(i,j) - &
        dt*gamma_oc*kc*cTOC(i, j) - dt*kb/A*dAdy + &
        dt*B/A*(u_in*DO_in - u_out*cDO(i, j))- V*v1DO + v2DO

        if (i == N) then
            cDO(i, j+1)= cDO(i, j) + dt*gamma_oy*(Gy - Ry)*cChla(i,j) - &
            dt*gamma_oc*kc*cTOC(i, j) - dt*kb/A*dAdy + &
            dt*B/A*(u_in*do_in - u_out*cDO(i, j))+ aeration(Vw, cDO(i,j), dz, t(i,j)) + 2*s2*(cDO(i-1,j) -cDO(i, j))
        end if 
        
        if (cDO(i, j+1) < 0 ) then
            cDO(i, j+1) = 0
        endif

        ! SS
        v1SS = s1/2*(cSS(i+1, j)- cSS(i-1, j))
        v2SS = s2*(cSS(i+1, j) - 2*cSS(i,j) + cSS(i-1, j))

        cSS(i, j+1) = cSS(i, j) + W0*v1ss + &
        dt*B/A*(u_in*ss_in - u_out*cSS(i, j))- V*v1SS + v2SS

        if (i == N) then
            cSS(i, j+1) = cSS(i, j) + W0*v1ss + &
            dt*B/A*(u_in*ss_in - u_out*cSS(i, j)) + 2*s2*(cSS(i-1,j) - cSS(i, j)) 
        end if

        if (cSS(i, j+1) < 0 ) then
            cSS(i, j+1) = 0
        endif
        
        ! Temperature
        v1T = s1/2*(t(i+1, j)- t(i-1, j))
        v2T = s2*(t(i+1, j) - 2*t(i,j) + t(i-1, j))

        t(i, j+1) = t(i, j) + dt*B/A*(u_in*T_in - u_out*t(i, j))- V*v1T + v2T - &
        dt/(rho*C_heat)*yeta*phi_y

        if (i == N) then
            t(i, j+1) = t(i, j) + dt*B/A*(u_in*T_in - u_out*t(i, j)) + 2*s2*(t(i-1, j) - t(i, j)) - &
            dt/(rho*C_heat)*yeta*phi_y  !deviation of phi: yeta*phi_y
        end if

        if (t(i, j+1) < 0 ) then
            t(i, j+1) = 0
        endif

    end do

    open(10,file='output_chla.csv')
    open(11,file='output_toc.csv')
    open(12,file='output_tp.csv')
    open(13,file='output_tn.csv')
    open(14,file='output_do.csv')
    open(15,file='output_ss.csv')
    open(16,file='output_t.csv')

    if (j == 0) then
        write(10, 100) j*dt, cChla(0:N, 0)
        write(11, 100) j*dt, cTOC(0:N, 0)
        write(12, 100) j*dt, cTP(0:N, 0)
        write(13, 100) j*dt, cTN(0:N, 0)
        write(14, 100) j*dt, cDO(0:N, 0)
        write(15, 100) j*dt, cSS(0:N, 0)
        write(16, 100) j*dt, t(0:N, 0)
    else if (mod(j, 5) .eq. 0) then
        write(10, 100) j*dt, cChla(0:N, j)
        write(11, 100) j*dt, cTOC(0:N, j)
        write(12, 100) j*dt, cTP(0:N, j)
        write(13, 100) j*dt, cTN(0:N, j)
        write(14, 100) j*dt, cDO(0:N, j)
        write(15, 100) j*dt, cSS(0:N, j)
        write(16, 100) j*dt, t(0:N, j)
    end if

    end do

100 format (41(F13.6, ','), F13.6)
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    contains    
    
    real function light_intensity (elevation, phi_n)
    real, intent (in) :: elevation, phi_n

    real :: beta = 0.5, &                  ! Absorption coefficient [-]
            alpha_r = 0.06, &              ! Reflection coefficient [-]
            yeta = 0.4, &                  ! Light reduction coefficient [1/m]
            L = 80.0                       ! Depth [m]

        light_intensity = (1-beta)*(1-alpha_r)*phi_n*exp(-yeta*(L-elevation))
        
    end function light_intensity
    

    real function aeration (Vw, cDO, dz, wt)
    real, intent (in) :: Vw, cDO, dz, wt

    real :: Vw0 = 2.0, &                     ! Wind velocity [m/s]
            Lx0 = 440.0, &                   ! Film thickness [um]
            Lx, DO_sat, x, mdo

        if (Vw >= Vw0) then
            Lx = Lx0*(Vw/Vw0)**(-2)
        else 
            Lx = Lx0
        end if 

        DO_sat = 14.652-0.41022*wt+0.007991*wt**2-0.000077774*wt**3
        x = 0.024716*(wt - 20.0)
        mdo = 2*EXP(x)/1000000000            ! molecular diffusion coefficients [m2/s] 

        aeration = mdo/(Lx*dz)*(DO_sat - cDO)

    end function aeration


end program lake
