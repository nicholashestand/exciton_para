OBJDIR=./obj
SRCDIR=./src
g09SRCFILES=g09_calcInfo.f95 g09_hessian.f95 g09_mocoeff.f95 \
			g09_overlap.f95 g09_cicoeff.f95
g09OBJFILES=${g09SRCFILES:.f95=.o}
DEP=g09_commonvar.mod
HRSRCFILE=g09_lambda.f95
HROBJFILE=g09_lambda.o
HREXE=g09_lambda.exe
CTSRCFILE=g09_ctint.f95
CTOBJFILE=g09_ctint.o
CTEXE=g09_ctint.exe
TQSRCFILE=g09_tq.f95
TQOBJFILE=g09_tq.o
TQEXE=g09_tq.exe
vpath %.o ${OBJDIR}
vpath %.mod ${OBJDIR}
vpath %.f95 ${SRCDIR}
FLINKER = gfortran
LIBLINKER = -lblas -llapack

all : ${CTEXE} ${HREXE} ${TQEXE}

${TQEXE} : ${TQOBJFILE} ${g09OBJFILES}
	${FLINKER} $(addprefix ${OBJDIR}/,$(^F)) ${OBJDIR}/${DEP:.mod=.o} \
		-o $@ -I${OBJDIR} ${LIBLINKER}

${CTEXE} : ${CTOBJFILE} ${g09OBJFILES}
	${FLINKER} $(addprefix ${OBJDIR}/,$(^F)) ${OBJDIR}/${DEP:.mod=.o} \
		-o $@ -I${OBJDIR} ${LIBLINKER}

${HREXE} : ${HROBJFILE} ${g09OBJFILES}
	${FLINKER} $(addprefix ${OBJDIR}/,$(^F)) ${OBJDIR}/${DEP:.mod=.o} \
		-o $@ -I${OBJDIR} ${LIBLINKER}

%.o : %.f95 ${DEP}
	${FLINKER} -c ${SRCDIR}/$(<F) -o ${OBJDIR}/$(@F) -I${OBJDIR}

${DEP} : ${DEP:.mod=.f95}
	${FLINKER} -c $< -o ${OBJDIR}/${@:.mod=.o} -J${OBJDIR}

clean :
	@echo cleaning up
	@rm ${OBJDIR}/*.o ${OBJDIR}/*.mod 2>/dev/null || true
