!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
! Copyright (C) 2021-2023 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov
!                         Dhruv Desai, Sergio Pineda Flores, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!
! Maintenance:
!===============================================================================

module pert_const
   use kinds, only: dp
   implicit none
   public
   !double precision
   !integer,  parameter :: dp = kind(1.0d0)
   !to be consistent with QE, set dp to be the same as QE/kinds
   !integer,  parameter :: dp = selected_real_kind(14,200)
   !constants
   real(dp), parameter :: pi = 3.141592653589793238462643383279_dp
   real(dp), parameter :: twopi = 2.0_dp*pi
   real(dp), parameter :: ryd2ev  = 13.605698066_dp
   real(dp), parameter :: ryd2mev = ryd2ev*1.0E3_dp
   ! Wu Li's paper also uses this convention(mev2hz_hbar), toTHz = 20670.687 ~  
   ! if you want Thz, please use mev2thz_h, 
   ! if you want rad/picosecond, please use mev2hz_hbar
   real(dp), parameter :: mev2hz_hbar  = 1.519267296E12_dp
   real(dp), parameter :: ryd2hz_hbar   = ryd2ev * mev2hz_hbar * 1E+3
   real(dp), parameter :: ryd2thz_hbar  = ryd2hz_hbar * 1E-12 ! ryd2thz_hbar = toTHz = 20670.687       
   real(dp), parameter :: mev2hz_h  = 2.41798924262E11_dp 
   real(dp), parameter :: ryd2hz_h   = ryd2ev * mev2hz_h * 1E+3
   real(dp), parameter :: ryd2thz_h  = ryd2hz_h * 1E-12       
   real(dp), parameter :: kelvin2eV = 8.6173427909E-05_dp
   real(dp), parameter :: bohr2ang = 0.52917721092_dp
   real(dp), parameter :: e2 = 2.0_dp  !e^2 in Rydberg atomic unit
   real(dp), parameter :: timeunit = 4.8377687E-17_dp !Rydberg atomic unit t0 in s
   real(dp), parameter :: unitcharge = 1.60217733E-19_dp !electron charge e in C, if so the q^2 could be treat as 1. Otherwise if we define unitcharge = 1.1329105*10^-19, then q^2 = 2.0
   real(dp), parameter :: tesla2ryd = 4.254381162420723E-6_dp !Convert tesla to ryd
   real(dp), parameter :: efield_tolerance = 1.0E-12_dp ! tolerance for the absolute value of external electric field
   real(dp), parameter :: drift_accel_tolerance = 1.0E-12_dp ! tolerance for the drift acceleration in a.u.
   
   complex(dp), parameter :: czero = (0.0_dp,0.0_dp)
   complex(dp), parameter :: cone  = (1.0_dp,0.0_dp)
   complex(dp), parameter :: ci    = (0.0_dp,1.0_dp)
end module pert_const
