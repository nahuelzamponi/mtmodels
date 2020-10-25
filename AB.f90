include 'random_module.f90'

program mito3
!
! Programa de generacion de red mitocondrial
! Noviembre de 2016
! version 1 - Calcula distribucion de tamaños de clusters

USE randomize
USE omp_lib

!----------------------------------------------------------------
implicit none

integer,parameter                   :: L = 10000 ! número de dimeros

!------------   parametros fijos   ------------------------------

integer,parameter                   :: N = 2*L
real*8,parameter                    :: PI = 3.14159265359

!------------ variables -----------------------------------------

integer					:: i, j, k, ii, jj, kk, idum, mmm
integer					:: iteraciones, muestras
integer					:: num_clusters, max_cluster, max_cluster_size
integer, allocatable, dimension(:,:)	:: equiv, vec
integer*1, allocatable, dimension(:)	:: grado, rgrado
integer, allocatable, dimension(:)	:: par, rpar, N1, N2, N3, nuevo, cluster, list
integer*4, allocatable, dimension(:)	:: size_cluster
real*4, allocatable, dimension(:)	:: Ns
real*4					:: Na
integer					:: N1max, N2max, N3max, Nmax
real					:: c1, c2
real					:: grado_m, giant_m, smedio, smedio_m
character(40)				:: file1, file2
character(20)				:: str

!---------------- archivos -------------------------------
! file1='hL10000c10p01-c20p000009.dat'
file1 = 'salida_prueba.dat'
!
!---------------- parametros de simulacion ---------------
!
idum = 28373466      ! semilla generador de numeros aleatorios
iteraciones = 3*L
c1 = 0.01
c2 = 0.000009
muestras = 1000
!---------------- INICIALIZACION --------------------------

mmm = mzranset(521288629, 362436069, 16163801, idum)

allocate(vec(N,3))
allocate(equiv(N,2))
allocate(grado(N))
allocate(rgrado(N))
allocate(par(N))
allocate(nuevo(N))
allocate(list(N))
allocate(cluster(N))
allocate(size_cluster(N))
allocate(N1(N))
allocate(N2(N))
allocate(N3(N))
allocate(Ns(N))

!---------------- PROGRAMA PRINCIPAL ------------------

Ns = 0
do k = 1, muestras
	call inicializa_arrays

	do ii=1, iteraciones  !transitorio
		call reaccion
	end do

	print *, k	

	call genera_red
	call get_clusters
	do jj = 1, num_clusters
		Ns(size_cluster(jj)) = Ns(size_cluster(jj)) + 1/float(Nmax)
	end do
end do

open(unit=1, file=file1)
do ii=1, max_cluster_size
	!if (Ns(ii).NE.0) then
		!write(1,*) ii,' ',Ns(ii)/float(muestras)
	!end if
	Na=0
	do jj=ii,max_cluster_size
		Na = Na + Ns(jj)/float(muestras) !cumulative distribution
	end do
	write(1, *) ii, ' ', Na
end do
close(1)

!------------------------------------------------------

deallocate(vec)
deallocate(grado)
deallocate(equiv)
deallocate(rgrado)
deallocate(par)
deallocate(nuevo)
deallocate(list)
deallocate(cluster)
deallocate(size_cluster)
deallocate(N1)
deallocate(N2)
deallocate(N3)
deallocate(Ns)

!------------------------------------------------------

contains

subroutine inicializa_arrays

integer	i

equiv=0

do i=1,N
	grado(i)=1
end do

do i=1,N-1,2
	par(i)=i+1
end do

do i=2,N,2
	par(i)=i-1
end do   

do i=1,N
	N1(i)=i
end do 

N2 = 0
N3 = 0

N1max = N
N2max = 0
N3max = 0
Nmax = 0
num_clusters = 0

end subroutine inicializa_arrays

!------------------------------------------------------

subroutine reaccion

real	z, zt, alfa1, alfa2, alfa3, alfa4, alfa12, alfa123, Rt

alfa1 = c1*N1max*(N1max - 1)/2.
alfa2 = float(N2max)
alfa3 = c2*N1max*N2max
alfa4 = 3*N3max/2.
alfa12 = alfa1 + alfa2
alfa123 = alfa12 + alfa3
Rt = alfa123 + alfa4

z=random()
zt= z*Rt

if (zt > alfa123) then !fision tip to side
	call reaccion4
else if (zt > alfa12) then !fusion tip to side
   call reaccion3
else if (zt > alfa1) then !fision tip to tip
	call reaccion2
else !fusion tip to tip
	call reaccion1
end if

end subroutine reaccion

!------------------------------------------------------------------------

subroutine reaccion1 !fusion tip to tip

integer	ii, indice1, indice2, sitio1, sitio2
logical	salida

if (N1max > 2) then
	indice1 = int(random()*N1max + 1)
	sitio1 = N1(indice1)
	do ii = indice1, N1max - 1
		N1(ii) = N1(ii + 1)
	end do
	N1max = N1max - 1
	if (grado(par(sitio1)) == 1) then
		salida = .true.
		do while (salida)
			indice2 = int(random()*N1max + 1)
			if (N1(indice2) /= par(indice1)) then
				salida = .false.
			end if 
		end do
	else
	indice2 = int(random()*N1max + 1)
	end if

	sitio2 = N1(indice2)
	do ii = indice2, N1max - 1
		N1(ii) = N1(ii + 1)
	end do
	N1max = N1max - 1
	N2max = N2max + 1
	N2(N2max) = sitio1
	equiv(sitio1, 1) = sitio2
	equiv(sitio2, 1) = sitio1
	grado(sitio1) = 2
	grado(sitio2) = 2
end if

end subroutine reaccion1

!------------------------------------------------------------------------

subroutine reaccion2 !fision tip to tip

integer  ii, jj, indice1, indice2, sitio1, sitio2
logical  salida

if (N2max >= 1) then
	indice1 = int(random()*N2max + 1)
	sitio1 = N2(indice1)
	do ii = indice1, N2max - 1
		N2(ii) = N2(ii + 1)
	end do
	N2max = N2max - 1
	sitio2 = equiv(sitio1, 1)
	N1max = N1max+2
	N1(N1max - 1) = sitio1
	N1(N1max) = sitio2
	equiv(sitio1, 1) = 0
	equiv(sitio2, 1) = 0
	grado(sitio1) = 1
	grado(sitio2) = 1
end if

end subroutine reaccion2

!------------------------------------------------------------------------

subroutine reaccion3 !fusion tip to side

integer  ii,indice1,indice2,indice3,sitio1,sitio2,sitio3

if ((N1max >= 1).and.(N2max >= 1)) then
	indice1 = int(random()*N1max + 1)
	sitio1 = N1(indice1)
	indice2 = int(random()*N2max + 1)
	sitio2 = N2(indice2)
	sitio3 = equiv(sitio2,1)
	do ii = indice1, N1max - 1
		N1(ii) = N1(ii + 1)
	end do
	N1max = N1max - 1
	do ii = indice2, N2max - 1
		N2(ii) = N2(ii + 1)
	end do
	N2max = N2max - 1
	N3max = N3max + 1
	N3(N3max) = sitio1
	equiv(sitio1, 1) = sitio2
	equiv(sitio1, 2) = sitio3
	equiv(sitio2, 1) = sitio1
	equiv(sitio2, 2) = sitio3
	equiv(sitio3, 1) = sitio1
	equiv(sitio3, 2) = sitio2
	grado(sitio1) = 3
	grado(sitio2) = 3
	grado(sitio3) = 3
end if

end subroutine reaccion3

!------------------------------------------------------------------------

subroutine reaccion4 !fision tip to side

integer  ii, indice1, sitio1, sitio2, sitio3

if (N3max >= 1) then
	indice1 = int(random()*N3max + 1)
	sitio1 = N3(indice1)
	sitio2 = equiv(sitio1, 1)
	sitio3 = equiv(sitio1, 2)
	do ii = indice1, N3max - 1
		N3(ii) = N3(ii + 1)
	end do
	N3max = N3max - 1
	N1max = N1max + 1
	N1(N1max) = sitio1
	N2max = N2max + 1
	N2(N2max) = sitio2
	equiv(sitio1, 1) = 0
	equiv(sitio1, 2) = 0
	equiv(sitio2, 1) = sitio3
	equiv(sitio2, 2) = 0
	equiv(sitio3, 1) = sitio2
	equiv(sitio3, 2) = 0
	grado(sitio1) = 1
	grado(sitio2) = 2
	grado(sitio3) = 2
end if

end subroutine reaccion4

!------------------------------------------------------------------------
subroutine genera_red

integer ii, jj

vec = 0
rgrado = 0

!construye red de vecinos
nuevo = 0
jj = 0
do ii = 1, N1max
	jj = jj + 1
	nuevo(N1(ii)) = jj
	rgrado(jj) = 1
end do

do ii = 1, N2max
	jj = jj + 1
	nuevo(N2(ii)) = jj
	nuevo(equiv(N2(ii), 1)) = jj
	rgrado(jj) = 2
end do

do ii = 1, N3max
	jj = jj + 1
	nuevo(N3(ii)) = jj
	nuevo(equiv(N3(ii),1)) = jj
	nuevo(equiv(N3(ii),2)) = jj
	rgrado(jj) = 3
end do
Nmax = jj

jj = 0
do ii = 1, N1max
	jj = jj + 1
	vec(jj, 1) = nuevo(par(N1(ii)))
end do

do ii = 1, N2max
	jj = jj + 1
	vec(jj, 1) = nuevo(par(N2(ii)))
	vec(jj, 2) = nuevo(par(equiv(N2(ii),1)))
end do

do ii = 1, N3max
	jj = jj + 1
	vec(jj, 1) = nuevo(par(N3(ii)))
	vec(jj, 2) = nuevo(par(equiv(N3(ii), 1)))
	vec(jj, 3) = nuevo(par(equiv(N3(ii), 2)))
end do

end subroutine genera_red
!------------------------------------------------------------------------

subroutine graba_red(ii) !graba red

integer ii, jj, kk

write(str, '(I4)') ii
file1 = 'red-'//trim(adjustl(str))//'.dat'
open(unit = 2, file = file1)
do jj = 1, Nmax
	do kk = 1, rgrado(jj)
	write(2, *)jj, '  ', vec(jj, kk)
	end do
end do   
close(1)

end subroutine graba_red

!------------------------------------------------------------------------

subroutine graba_pajek(ii) !graba red  en formato pajek

integer ii, jj, kk

write(str, '(I4)') ii
file1 = 'red-'//trim(adjustl(str))//'.net'
open(unit = 2, file = file1)
write(2,*)'*Vertices ', Nmax
do jj = 1, Nmax
	write(2,*)jj, ' ', '"', jj, '"'
end do
write(2, *)'*Arcs'
do jj = 1, Nmax
	do kk = 1, rgrado(jj)
		write(2, *)jj, '  ', vec(jj, kk), ' 1'
	end do
end do   
close(1)

end subroutine graba_pajek

!------------------------------------------------------------------------

subroutine genera_cluster(ii, color)
! genera cluster a partir del sitio ii

integer color, i, ii, maxlist, sitio

list = 0
sitio = ii
cluster(sitio) = color
maxlist = 1
list(maxlist) = sitio
i = 1
do while (i <= maxlist)
	if (rgrado(list(i)) == 1) then
		sitio = vec(list(i), 1)
		if(cluster(sitio) == 0) then
			cluster(sitio) = color
			maxlist = maxlist + 1
			list(maxlist) = sitio
		end if
	else if (rgrado(list(i)) == 2) then
		sitio = vec(list(i), 1)
		if(cluster(sitio) == 0) then
			cluster(sitio) = color
			maxlist = maxlist + 1
			list(maxlist) = sitio
		end if
		sitio = vec(list(i), 2)
		if(cluster(sitio) == 0) then
			cluster(sitio) = color
			maxlist = maxlist + 1
			list(maxlist) = sitio
		end if
	else
		sitio = vec(list(i), 1)
		if(cluster(sitio) == 0) then
			cluster(sitio) = color
			maxlist = maxlist + 1
			list(maxlist) = sitio
		end if
		sitio = vec(list(i), 2)
		if(cluster(sitio) == 0) then
			cluster(sitio) = color
			maxlist = maxlist + 1
			list(maxlist) = sitio
		end if
		sitio = vec(list(i), 3)
		if(cluster(sitio) == 0) then
			cluster(sitio) = color
			maxlist = maxlist + 1
			list(maxlist) = sitio
		end if
	end if
  	i = i + 1
end do
end subroutine genera_cluster
!------------------------------------------------------------------------

subroutine get_clusters

integer color, ii, jj

cluster = 0
color = 0

do ii = 1, Nmax
	if (cluster(ii) == 0) then
		color = color + 1
		call genera_cluster(ii, color)
	end if
end do

num_clusters = color
size_cluster = 0
 
do jj = 1, num_clusters !extrae clusters de la red
	do ii = 1, Nmax
		if (cluster(ii) == jj) then
			size_cluster(jj) = size_cluster(jj) + 1
		end if
	end do
end do

max_cluster_size=size_cluster(1) !extrae el cluster gigante
max_cluster = 1
do ii = 2, num_clusters
	if (size_cluster(ii) > max_cluster_size) then
		max_cluster_size=size_cluster(ii)
		max_cluster = ii
	end if
end do

end subroutine get_clusters

!------------------------------------------------------------------------

subroutine average_clusters_size

integer ii
real    sum, sum2

do ii = max_cluster, num_clusters - 1 !extrae cluster gigante
	size_cluster(ii) = size_cluster(ii + 1)
end do
num_clusters = num_clusters - 1

sum = 0
sum2 = 0
do ii = 1, num_clusters
	sum = sum + size_cluster(ii)/float(Nmax)
	sum2 = sum2 + size_cluster(ii)**2/float(Nmax)
end do

smedio = sum2/sum

end subroutine average_clusters_size

!---------------------

end program
