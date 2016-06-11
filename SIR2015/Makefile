# ****************************************************************************
# 
# Fichero   : Makefile
# Proyecto  : Metodos de inversion
# Dpto.     : Departamento de Astrofisica - ULL
# Autor     : Basilio Ruiz Cobo - brc@ll.iac.es
# Proposito : Mantenimiento de los ejecutables de inversion
#
# Uso       : > make [fc=compilador]               ... genera sir.x
#             > make [fc=compilador] lande.x       ... genera lande.x
#             > make [fc=compilador] modelador2.x  ... genera modelador2.x
#
#             donde fc puede ser gfortran, ifort o nada (defecto = ifort).
#             Por ejemplo:
#
#             > make sir.x
#             > make fc=gfortran lande.x
#             > make fc=ifort modelador2.x
#
# 27/01/04  : primera version.
# 12/11/15  : el compilador a usar puede entrar por linea de comandos.
#
# ****************************************************************************


# Recordatorio de la ultima compilacion.

.KEEP_STATE:

# Compilador, flags del compilador y flags de construccion de la libreria.
#
# Extension de lineas de codigo: -ffixed-line-length-none o -extend_source

ifeq ($(fc),gfortran)
FC = gfortran
FFLAGS = -O3 -ffixed-line-length-none -fno-automatic -fbounds-check -ffast-math -funroll-loops
else ifeq ($(fc),ifort)
FC = ifort
FFLAGS = -O2 -extend_source
else
FC = ifort
FFLAGS = -O2 -extend_source
endif
ARFLAGS = ruv
 
# Fichero de configuracion (definicion de dimensiones).

FILECONF = PARAMETER
FILEMAKE = .make.lib

# Ficheros fuentes de los ejecutables.

MAIN1SRC = sir.f
MAIN2SRC = modelador2.f
MAIN3SRC = lande.f
MAIN4SRC = modelador3.f

# Libreria de inversion y todos sus modulos.

LIB = libreria.a

LIBSRC = alphasub.f \
	amp2.f \
	atmdatb.f \
	automatico.f \
	automaticop.f \
	blends2.f \
	blendscon2.f \
	cabecera.f \
	chachi.f \
	comprime2.f \
	conhsra.f \
	deconv.f \
	deconv2.f \
	delta1.f \
	deriva.f \
	deriva4cuad.f \
	derivacuad.f \
	dsaha.f \
	dtplanck.f \
	equi.f \
	equisubmu.f \
	equisubmu_cont.f \
	equisubmu_contmag.f \
	escribeFR.f \
	fint.f \
	fperfil2.f \
	fperfil2err.f \
	gasb.f \
	gase.f \
	hermite.f \
	hermite_c.f \
	inicia_pefrompgt.f\
	intmodel.f \
	kappach.f \
	leeabun.f \
	leelineasii.f \
	leemalla2.f \
	leemallab.f \
	leemodi22.f \
	leemodi222.f \
	leeuve2.f \
	leeuveobs.f \
	leeuveobsindic.f \
	locate.f \
	lubksb.f \
	ludcmp.f \
	marqcoef2.f \
	marqcoeferr.f \
	marquarderr.f \
	marquardt2.f \
	matinx.f \
	matinx2.f \
	mensaje.f \
	meves.f \
	model.f \
	molecb.f \
	mreadi2.f \
	mreadi3.f \
	mreadi33.f \
	mvoigt.f \
	mvoigtc.f \
	nodos2.f \
	nodos2aut.f \
	penalty.f \
	penalty2.f \
    pgpefromrho.f \
	petaufrompgzro.f \
	pemufrompgtaut.f \
	pgzrofrompetau.f \
	planck.f \
	plank2.f \
	polint.f \
	reduce2.f \
	refrax.f \
	rnorma.f \
	saha.f \
	splinb.f \
	splines.f \
	splines22.f \
	sublande.f \
	sublande2.f \
	svbksb.f \
	svdcmp.f \
	svdcmp2.f \
	svdmatriz2.f \
	taulinea.f \
	taulinea2.f \
	thermosub.f \
	gasb_thermo.f \
	voigt.f \
	zeeman.f \
	zeeman_jk.f

# Referencias para el reemplazamiento de sufijos.

MAIN1OBJ = ${MAIN1SRC:.f=.o}
MAIN2OBJ = ${MAIN2SRC:.f=.o}
MAIN3OBJ = ${MAIN3SRC:.f=.o}
MAIN4OBJ = ${MAIN4SRC:.f=.o}

LIBOBJ = ${LIBSRC:.f=.o}


#
# A)
# Ejecutable sir.x ('uso> make' o 'uso> make sir.x').
#

sir.x: ${MAIN1OBJ} ${LIB} ${FILEMAKE}
	${LINK.f} -o $@ ${MAIN1OBJ} ${LIB}


#
# B)
# Ejecutable modelador2.x ('uso> make modelador2.x').
#

modelador2.x: ${MAIN2OBJ} ${LIB}
	${LINK.f} -o $@ ${MAIN2OBJ} ${LIB}

#
# C)
# Ejecutable lande.x ('uso> make lande.x').
#

lande.x: ${MAIN3OBJ} ${LIB}
	${LINK.f} -o $@ ${MAIN3OBJ} ${LIB}

# B2)
# Ejecutable modelador3.x ('uso> make modelador3.x').
#

modelador3.x: ${MAIN4OBJ} ${LIB}
	${LINK.f} -o $@ ${MAIN4OBJ} ${LIB} 


#
# D)
# Actualizacion de la libreria.
#

${LIB}: ${LIBOBJ}
	${AR} ${ARFLAGS} $@ $?


#
# E)
# Reconstruccion de la libreria.
#

${FILEMAKE}: ${FILECONF}
	${COMPILE.f} ${LIBSRC}
	${AR} ${ARFLAGS} ${LIB} ${LIBOBJ}
	@touch $@


# ****************************************************************************