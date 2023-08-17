"""
      Calculation for control position -10.00
   Mass=       0.185 kg

   ___Center of Gravity Position - Body axis____
    CoG_x=      0.0600 m
    CoG_y=      0.0000 m
    CoG_z=      0.0000 m

   ___Inertia - Body Axis - CoG Origin____
    Ibxx=       0.004 kg.m²
    Ibyy=       0.004 kg.m²
    Ibzz=       0.006 kg.m²
    Ibxz=           0 kg.m²

         Rotating the elevator by -1.00°, total angle is -1.00°
      Creating the unit RHS vectors...
      Creating the influence matrix...
      Performing LU Matrix decomposition...
      Solving the LU system...
      Time for linear system solve: 0.080 s
      Searching for zero-moment angle... Alpha=1.16357°
      Creating source strengths...
      Calculating doublet strength...
      Calculating speed to balance the weight...VInf = 8.57402 m/s

      ___Inertia - Stability Axis - CoG Origin____
      Isxx=    0.004001 
      Isyy=       0.004 
      Iszz=    0.005999 
      Isxz=   -4.06e-05 

      Calculating the stability derivatives
         Creating the RHS translation vectors
         LU solving for RHS - longitudinal
         Calculating forces and derivatives - lateral
         Creating the RHS rotation vectors
         LU solving for RHS - lateral
         Calculating forces and derivatives - lateral
      Calculating the control derivatives

      Longitudinal derivatives
      Xu=   -0.028858         Cxu=   -0.039385
      Xw=      0.1091         Cxa=     0.14889
      Zu=    -0.42128         Czu=   0.0028013
      Zw=     -3.1536         CLa=      4.3039
      Zq=    -0.44812         CLq=      6.9449
      Mu= -0.00050102         Cmu=  -0.0038824
      Mw=   -0.081849         Cma=    -0.63425
      Mq=   -0.098813         Cmq=     -8.6951
      Neutral Point position=   0.08595 m


      Lateral derivatives
      Yv=    -0.61866         CYb=    -0.84433
      Yp=    0.011975         CYp=    0.040857
      Yr=     0.26989         CYr=     0.92085
      Lv=   -0.044646         Clb=   -0.076165
      Lp=   -0.090387         Clp=    -0.38549
      Lr=    0.019467         Clr=    0.083026
      Nv=     0.21312         Cnb=     0.36358
      Np=   -0.037193         Cnp=    -0.15862
      Nr=   -0.097294         Cnr=    -0.41495

      Control derivatives 
      Xde=   -0.077213        CXde=    -0.01229
      Yde= -4.8915e-06        CYde=  -7.786e-07
      Zde=     -5.8507        CZde=    -0.93128
      Lde= -2.4741e-08        CLde= -4.9226e-09
      Mde=     -1.9823        CMde=     -1.7916
      Nde= -3.8288e-06        CNde= -7.6181e-07


      _____State matrices__________
       Longitudinal state matrix
             -0.155991             0.58972                   0               -9.81
              -2.27721            -17.0465             6.15176                   0
             -0.125254            -20.4623            -24.7031                   0
                     0                   0                   1                   0
       Lateral state matrix
              -3.34412           0.0647292            -7.11515                9.81
              -11.5206            -22.5308             5.03081                   0
               35.6033             -6.0472             -16.252                   0
                     0                   1                   0                   0

      _____Control Matrices__________
       Longitudinal control matrix
         -0.4173653
          -31.62547
          -495.5809
                  0

       Lateral control matrix
      -2.644068e-05
       2.934983e-07
      -0.0006382265
                  0



      ___Longitudinal modes____

      Eigenvalue:     -20.88+   -10.57i   |      -20.88+    10.57i   |    -0.07713+  -0.8891i   |    -0.07713+   0.8891i
                    _____________________________________________________________________________________________________
      Eigenvector:         1+        0i   |           1+        0i   |           1+        0i   |           1+        0i
                       17.58+   -40.79i   |       17.58+    40.79i   |     -0.1045+-0.003387i   |     -0.1045+ 0.003387i
                      -80.68+   -4.821i   |      -80.68+    4.821i   |     0.08151+ 0.005757i   |     0.08151+-0.005757i
                       3.169+   -1.374i   |       3.169+    1.374i   |    -0.01432+  0.09043i   |    -0.01432+ -0.09043i



      ___Lateral modes____

      Eigenvalue:     -22.54+        0i   |      -9.786+   -15.59i   |      -9.786+    15.59i   |    -0.01043+        0i
                    _____________________________________________________________________________________________________
      Eigenvector:         1+        0i   |           1+        0i   |           1+        0i   |           1+        0i
                       8.248+        0i   |     -0.6322+  0.07478i   |     -0.6322+ -0.07478i   |    -0.02019+        0i
                       2.269+        0i   |      0.9201+    2.148i   |      0.9201+   -2.148i   |         2.2+        0i
                     -0.3659+        0i   |     0.01483+ -0.03125i   |     0.01483+  0.03125i   |       1.935+        0i

      Calculating aerodynamic coefficients in the far field plane
        Calculating point    1.16°....
      Computing On-Body Speeds...
      Computing Plane for alpha=   1.16°
       Calculating aerodynamic coefficients...
         Calculating wing...Main Wing
         Calculating wing...Elevator
         Calculating wing...Fin

   Phillips formulae:
       Phugoid eigenvalue:      -0.02686+  0.87538i
               frequency:  0.139 Hz
               damping:    0.031
       Dutch-Roll eigenvalue:   -9.75549+ 15.57029i
               frequency:  2.924 Hz
               damping:    0.627
"""



"""根軌跡を描画するプログラム"""
from control.matlab import *
from matplotlib import pyplot as plt
import numpy as np
plt.rc('font', size=12)

rho = 1.2
S = 0.19
Vh = 0.5
aw = 2*3.14*0.8
at = 2*3.14*0.6
lt = 0.4
U = 8.57402
c = 0.176
h = 0.33
hnw = 0.25
tau = 0.3
I = 0.0041

Ma = 0.5*rho*U*U*S*c*( aw*(hnw - h) + at*Vh )/I
Mde = 0.5*rho*U*U*S*c*at*Vh*tau/I
Mq = 0.5*rho*U*S*c*at*Vh*lt/I

omega_a = 2*3.14*3
zeta_a = 2*0.8*omega_a

Kd = 7/100
Ki = 0.0
K = 2.0

#Kd = 7/100
#Ki = 0.8
#K = 0.6

#計算ゲイン範囲
gain_range = np.linspace(0,K,1000)

T = [[1+0j,1+0j,1+0j,1+0j],
                [17.58-40.79j,17.58+40.79j,-0.1045-0.003387j,-0.1045+0.003387j],
                [-80.68-4.821j,-80.68+4.821j,0.08151+0.005757j,0.08151-0.005757j],
                [3.169-1.374j,3.169+1.374j,-0.01432+0.09043j,-0.01432-0.09043j]]
T = np.array(T)
LAMBDA = np.diag([-20.88-10.57j,-20.88+10.57j,-0.07713-0.8891j,-0.07713+0.8891j])

B =  np.array([[-0.4173653],
                 [-31.62547],
                 [-495.5809],
                 [0]])

                 

A = T@LAMBDA@np.linalg.inv(T)


#空力モデル
#sys = ss(A,B,np.array([0,0,0,1]),0)
#sys = tf(sys)
sys = tf([495.6,7880,1899],[1,41.91,554.9,117.7,436.2])
print(sys)
#サーボモデル
servo = tf([omega_a**2],[1,zeta_a,omega_a**2])
#制御モデル
controltf = tf([Kd,0],[1]) + tf([1],[1]) + tf([Ki],[1,0])
#システム全体
sysall = sys*servo*controltf
print(sys)
#根軌跡
roots_all, gains_all = rlocus(sysall,xlim=[-50,30],ylim=[-30,30],kvect=gain_range)
print(roots_all)
print(gains_all)
plt.savefig("long/rlocus_all.png")
plt.savefig("long/rlocus_all.eps")
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.savefig("long/rlocus_all_phu.png")
plt.savefig("long/rlocus_all_phu.eps")
plt.clf()

roots, gains = rlocus(sys*controltf,xlim=[-50,30],ylim=[-30,30],kvect=gain_range)
plt.savefig("long/rlocus.png")
plt.savefig("long/rlocus.eps")
print(roots)
print(gains)
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.savefig("long/rlocus_phu.png")
plt.savefig("long/rlocus_phu.eps")
plt.clf()

bode(sys,servo,controltf*K,sysall*K)
plt.savefig("long/bode.png")
plt.savefig("long/bode.eps")
plt.clf()

bode(sys,servo,controltf*K,sysall*K)
plt.savefig("long/bode_w.png")
plt.savefig("long/bode_w.eps")
plt.clf()

bode(sysall*K/(1+sysall*K),sys*controltf*K/(1+sys*controltf*K),initial_phase=0.0)
plt.savefig("long/bode_YR.png")
plt.savefig("long/bode_YR.eps")
plt.clf()