program panache_param
use netcdf
  implicit none

  integer :: ncid, netcdf_file,varid_b,varid_alt,varid_t, uid(2)
  integer, parameter :: nz = 100, nt = 1000
  real, parameter :: dt = 1.0d0, dz = 1.0d0, eps = 0.005d0
  real :: start_time, end_time, elapsed_time
  real :: b(nt, nz), btend_Q(nz), btend_param(nz)
  real :: Q, alpha, b_init(nz)
  integer :: it, iz
  real :: temps(nt), ZZ(nz)
  character(len=50) :: file_name

  ! Mesurer le temps avant une section de code
  ! call cpu_time(start_time)
  ! Initialisation des paramètres
  Q = 0.05
  alpha = 0.1
  file_name = "output_fortran.nc"

  ! Initialisation du tableau b
  ! Première partie de b_init : linéaire de 1 à 0.5
  do iz = 1, nz / 10
     b_init(iz) = 1.0d0 - 0.5d0 * (iz - 1) / (nz / 10 - 1)
  end do

  ! Seconde partie de b_init : linéaire de 0.5 à 2
  do iz = nz / 10 + 1, nz
     b_init(iz) = 0.5d0 + 1.5d0 * (iz - nz / 10 - 1) / (nz - nz / 10 - 1)
  end do

  ! Initialisation de b avec b_init
  b(1, :) = b_init

  ! Temps
  do it = 1, nt
     temps(it) = it * dt
  end do

  do iz = 1, nz
     ZZ(iz) = iz * dz
  end do

  ! Boucle temporelle
  do it = 1, nt - 1
     call Qf(nz, dt, Q, btend_Q)
     call Panache(b(it, :), nz, dt, dz, eps, alpha, btend_param)
     b(it+1, :) = b(it, :) + btend_Q + btend_param
     ! print *,b(it,50)
  end do

  ! Mesurer le temps après la section de code
  ! call cpu_time(end_time)

  ! Calculer le temps écoulé
  ! elapsed_time = end_time - start_time

  ! Afficher le temps écoulé
  ! print *, "Temps écoulé :", elapsed_time, "secondes"
  
  netcdf_file = nf90_create(file_name, NF90_CLOBBER,ncid)
  netcdf_file = nf90_def_dim(ncid, "Time", nt,uid(1))
  netcdf_file = nf90_def_dim(ncid, "Altitude",nz,uid(2))
  netcdf_file = nf90_def_var(ncid, "b", NF90_REAL,uid,varid_b)
  netcdf_file = nf90_def_var(ncid, "Altitude", NF90_REAL, uid(2), varid_alt)
  netcdf_file = nf90_def_var(ncid, "Time", NF90_REAL, uid(1), varid_t)
  netcdf_file = nf90_enddef(ncid)

  netcdf_file = nf90_put_var(ncid, varid_b, b)
  netcdf_file = nf90_put_var(ncid, varid_t, temps)
  netcdf_file = nf90_put_var(ncid, varid_alt, ZZ)
  netcdf_file = nf90_put_att(ncid, varid_t, "units", "s")
  netcdf_file = nf90_put_att(ncid, varid_alt, "units", "meters")
  netcdf_file = nf90_put_att(ncid, varid_b, "units", "ua")
  netcdf_file = nf90_put_att(ncid, varid_b, "long_name", "b")
  netcdf_file = nf90_close(ncid)

contains

  subroutine Qf(nz, dt, Q, Qf_out)
    integer, intent(in) :: nz
    real, intent(in) :: dt, Q
    real, intent(out) :: Qf_out(nz)
    integer :: iz

    do iz = 1, nz
       if (iz == 1) then
          Qf_out(iz) = Q * dt
       else
          Qf_out(iz) = -Q / (nz - 1) * dt
       end if
    end do
  end subroutine Qf

  subroutine Panache(b, nz, dt, dz, eps, alpha, Panache_out)
    real, intent(in) :: b(nz)
    integer, intent(in) :: nz
    real, intent(in) :: dt, dz, eps, alpha
    real, intent(out) :: Panache_out(nz)
    real :: f(nz), fb(nz), E(nz), D(nz), Gamma(nz), b_panache(nz)
    integer :: iz

    ! Initialisation
    f = 0.0d0
    fb = 0.0d0
    E = 0.0d0
    D = 0.0d0
    Gamma = 0.0d0
    b_panache = 0.0d0

    do iz = 1, nz
       if (iz == 1) then
          f(1) = 0.0d0
          if (nz > 1) then
             E(iz) = alpha / 2.0d0 * sqrt(max(0.0d0, b(iz) - b(iz+1))) * sqrt(dz)
          end if
          D(iz) = 0.0d0
          f(iz+1) = f(iz) + E(iz)
          b_panache(iz) = b(iz)
          fb(iz+1) = f(iz+1) * (b_panache(iz) - b(iz+1))
          Panache_out(iz) = -dt / dz * (fb(iz+1) - fb(iz))

       else if (iz < nz) then
          if (f(iz) > 0) then
             Gamma(iz) = alpha**2 * (b_panache(iz-1) - b(iz)) / (2.0d0 * f(iz))
          else
             Gamma(iz) = 0.0d0
          end if
          E(iz) = max(Gamma(iz), 0.0d0) + eps
          D(iz) = max(-Gamma(iz), 0.0d0) + eps
          f(iz+1) = max(f(iz) + E(iz) - D(iz), 0.0d0)
          if (f(iz+1) > 0) then
             b_panache(iz) = (f(iz) * b_panache(iz-1) + E(iz) * b(iz)) / (f(iz+1) + D(iz))
          else
             b_panache(iz) = b_panache(iz-1)
             E(iz) = 0.0d0
             D(iz) = f(iz)
          end if
          fb(iz+1) = f(iz+1) * (b_panache(iz) - b(iz+1))
          Panache_out(iz) = -dt / dz * (fb(iz+1) - fb(iz))

       else
          b_panache(iz) = b_panache(iz-1)
          E(iz) = 0.0d0
          D(iz) = f(iz)
          Panache_out(iz) = -dt / dz * (-fb(iz))
       end if
    end do
  end subroutine Panache

end program panache_param
