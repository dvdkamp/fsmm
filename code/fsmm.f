c     DV fuel moisture model foretran code


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialiser for parameter common block
      subroutine fmcpar(odeparms)
      external odeparms
      double precision parms(7)
      common /myparms/parms


      call odeparms(7, parms) 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialiser for forcing common block
      subroutine fmcforc(odeforcs)
      external odeforcs
      double precision forcs(8)
      common /myforcs/forcs
  
      call odeforcs(8, forcs)
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Rate of change and output variables
      subroutine fmccore (neq, t, y, ydot, out, IP)

      integer neq, IP(*)

      double precision t, y(neq), ydot(neq), out(*), 
     1     A,B,K_s,M_MAX,DR_O,r,L,
     2     L_out,L_up,Q_h,rh_s,psat_s,p_s,Ls,Q_e,
     3     m_all,t_all,Tsc,Fas,Fsc,P,Lv,c_wood, A_s,
     4     k_down_dir, k_down_diff, shape_factor, temp, 
     5     L_dn, rh, vap_press, precip, wind_speed,
     6     surf_area,surf_area_outer,volume,volume_o,volume_c,
     7     diff_factor,m_o,m_c,r_o_mid,r_c_mid,r_o,r_c,
     8     Re,Nu 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Constants
c
c     PI = 3.1415927 : a well-rounded number
c     RHO_S = 400   kg m-3 : Stick density  (taken from Nelson, 2000)
c     C_WATER = 4200 J kg-1 K-1 : water specific heat
c     RHO_A = 1.09266 kg/m^3 : Standard Density of Air. at 1175m ('76 Std. Atm.)
c     C_A = 1.00467e3 J kg-1 K-1  : Specific Heat of air From Stull (1988)
c     SB= 5.67e-8 W m-2 K-4 : Steph-Boltz constant
c     KELVIN = 273.15 K
c     M_W = 0.018 kg mol-1 : Molecular weight of water 
c     R_GAS = 8.314e-3 m^3 kPa K-1 mol-1 : Gas constant
c     M_FSP = 0.30 g/g: fibre saturation point (taken from Nelson, 2000)
c     E_G = 0.95: Forest floor emissivity
c     ALB_G = 0.185: Forest Floor albedo
c     ALB_S = 0.6, Stick Albedo
c     E_S = 0.85 Stick emissivity
c     G = 0.42 Specific Gravity (Simpson and TenWolde, 1999)

      REAL, PARAMETER :: PI = 3.1415927,RHO_S=400,
     1     C_WATER = 4200,RHO_A = 1.09266,C_A = 1.00467E3,SB= 5.67E-8,
     2     KELVIN = 273.15,M_W = 0.018, R_GAS = 8.314E-3,M_FSP = 0.3,
     3     E_G=0.95, ALB_G=0.165, ALB_S=0.6, E_S=0.85, G=0.42,
     4     k_a = 1.9E-5, v_a = 15.11E-6 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc BRING IN COMMON BLOCKS

c	Parameters
c    	A & B : constants for EMC equation
c	    K_s: Diffusivity (m^2 s^-1)
c	    M_MAX: Maximum allowable moisture content (%)
c	    DR_O: volume fraction of the stick taken up by the outer layer
c	    r: radius of stick (m)
c	    L: Length of stick (m)

c	    Forcing variables	 
c	    k_down_dir: direct solar radiation (W m^-2)
c	    k_down_diff: diffuse solar radiation (W m^-2)
c	    shape_factor: shape_factor for calculating captured direct solar radiation
c	    temp: temperature (C)
c	    rh: Relative Humidity (%)
c	    L_dn: downwelling longwave radiation (W m^-2)
c	    wind_speed: wind speed (m)
c	    precip: precipitation (mm)

      common /myparms/ A,B,K_s,M_MAX,DR_O,r,L
      common /myforcs/ k_down_dir, k_down_diff, shape_factor, temp, 
     1     L_dn, rh, wind_speed, precip


      if(IP(1) < 2) call rexit("nout should be at least 1")
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Dimensions

cc inner radius of outer layer
      r_o = r*SQRT(1 - DR_O)

cc middle radius of outer layer
      r_o_mid = (r_o + r)/2

cc radius of core 
      r_c = r_o

cc middle radius of core
      r_c_mid = r_c/2

      surf_area = 2*PI*r*L + 2*PI*r**2  
      surf_area_outer = 2*PI*r*L
      volume = PI*r**2*L 
      volume_o = volume - PI*r_o**2*L
      volume_c = PI*(r_c)**2*L
      diff_factor = PI*r*L + 2*PI*r**2
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Average moisture

      m_all = (y(2)  +  y(3))/(volume*RHO_S)

      t_all = (y(1)*volume_o + y(4)*volume_c)/volume

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc moistures: kg_h20 kg_wood^-1

      m_o = y(2)/(volume_o*RHO_S)
      m_c = y(3)/(volume_c*RHO_S)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc ENERGY BUDGET 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc EVAP AND SORPTION COOLING

      Re = wind_speed*r*2/v_a
      Nu = 0.17*Re**0.62
      res_h = r*2/(k_a*Nu)

c     surface relative humidity from Mathhews 2006 eq'n 16
      rh_s = exp(-4.19*M_W/(R_GAS*(y(1) + KELVIN))*exp(m_o*B + A))

c     ## saturation vapour pressure at stick surface temp
      if (y(1).gt.0) then
         psat_s = 0.611*exp(17.27*y(1)/(y(1) + 237.26))
      else   
         psat_s = 0.611*exp(21.87*y(1)/(y(1) + 265.6))
      endif

c     ## saturation vapour pressure at ambient temp
      if (temp.gt.0) then
         psat_a = 0.611*exp(17.27*temp/(temp + 237.26))
      else   
         psat_a = 0.611*exp(21.87*temp/(temp + 265.6))
      endif         
      
c      ## ambient vapour pressure
      vap_press = rh/100*psat_a
         
c     ## surface vapour pressure
      p_s = rh_s*psat_s
      
c     J / kg differential heat of sorptino
      Ls = 21000/M_W*exp(-14.*m_all)  

c     J/kg latent heat of evap (From Stull 1988)
      Lv = 2.501e6 - 0.00237e6*temp

c     ## Latent heat flux (W m-2) 
      Q_e = 1/res_h*(Lv+Ls)*M_W/(R_GAS*(temp+KELVIN))*(p_s-vap_press)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c LONGWAVE EMITTED (W m-2)

      L_out = E_S*SB*(y(1) + KELVIN)**4

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c upwelling longwave from ground (W m-2)

      L_up = E_G*SB*(temp + KELVIN)**4
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc sensible heat flux (W m-2)       
     
      Q_h = 1/res_h*(RHO_A*C_A)*(y(1) - temp)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Specific heat of moist wood (from Ragland et al., 1991)

      c_wood = 103.1 + 3.867*(y(1) + KELVIN)

      if (m_all < m_fsp) then
         A_s = (23.55*(y(1) + KELVIN) - 1320*m_all - 6191)*m_all
      else
         A_s = 0
      endif

      c_moist = (c_wood + m_all*C_WATER)/(1 + m_all) + A_s
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc thermal conductivity (W m^-2 K^-1) (Simpson and TenWolde, 1999, Eq'n 3-7)

       K_t = G*(0.1941 + 0.004064*m_all*100) + 0.01864

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Energy budget equation

c     Heat Conduction Flux into cylinder (in J s-1)
      Tsc = 2*PI*L*K_t*(y(1) - y(4))/log(r_o_mid/r_c_mid)

c     dTs/dt
      ydot(1) = 1/(c_moist*RHO_S*volume_o)* 
     1  (diff_factor*(E_S*(L_dn+L_up) +
     2               (1-ALB_S)*(k_down_diff+ALB_G*k_down_dir)) + 
     3  shape_factor*(1-ALB_S)*k_down_dir - 
     4  surf_area*(L_out + Q_h) - 
     5  surf_area_outer*Q_E - 
     6  Tsc)

      ydot(4) = 1/(c_moist*RHO_S*volume_c)*Tsc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Moisture Budgets

c     evap, kg m-2 s-1
      Fas = Q_e/(Lv + Ls)
  
c     diffusion to core: kg_h2o s-1
       Fsc = 2*PI*L*K_s*RHO_S*(m_o - m_c)/(log(r_o_mid/r_c_mid))

c     precipitation interception: kg_h20 hr-1
      P= min(precip*r*L,
     1     M_MAX*RHO_S*volume_o-
     2     max(0.0,(y(2)+(-Fas*surf_area_outer - Fsc)*3600)))

c     change in surface moisture, dm_s/dt : kg_h20 s-1 
       ydot(2) = -Fas*surf_area_outer- Fsc + P/3600

c     change in core moisture, dm_c/dt:  kg_h20 s-1
       ydot(3) = Fsc 

      out(1) = m_all
      out(2) = t_all

      return
      end


